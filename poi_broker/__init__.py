import logging
from datetime import datetime, timezone
from pathlib import Path

from flask import Flask, flash, jsonify, redirect, request, url_for
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from flask_login import LoginManager
from flask_sqlalchemy import SQLAlchemy
from flask_wtf.csrf import CSRFError, CSRFProtect
from sqlalchemy import event

from .extensions import huey
from .settings import build_app_config

#from werkzeug.middleware.profiler import ProfilerMiddleware
#import jinja2

# Initialize SQLAlchemy instance (outside create_app for import access)
db = SQLAlchemy()
login_manager = LoginManager()
csrf = CSRFProtect()
limiter = Limiter(key_func=get_remote_address, default_limits=[])

def _ensure_email_verification_expiry_column(app):
    """Add user.email_verification_token_expires on existing users DBs."""
    from sqlalchemy import inspect, text

    engine = db.engines.get('users', db.engine)
    inspector = inspect(engine)
    if 'user' not in inspector.get_table_names():
        return
    column_names = {col['name'] for col in inspector.get_columns('user')}
    if 'email_verification_token_expires' in column_names:
        return
    with engine.begin() as conn:
        conn.execute(text(
            'ALTER TABLE user ADD COLUMN email_verification_token_expires INTEGER'
        ))
    app.logger.info('Added user.email_verification_token_expires')


def _configure_sqlite_pragmas(dbapi_conn, connection_record):
    """Configure SQLite PRAGMAs for improved concurrent read performance."""
    cursor = dbapi_conn.cursor()
    cursor.execute("PRAGMA journal_mode = WAL") # Enable Write-Ahead Logging for better concurrency https://sqlite.org/wal.html
    cursor.execute("PRAGMA synchronous = NORMAL") # Balance durability and performance; WAL mode is safe with NORMAL https://sqlite.org/wal.html#durability
    cursor.execute("PRAGMA cache_size = -64000")  # 64MB in KB https://sqlite.org/pragma.html#pragma_cache_size
    cursor.execute("PRAGMA mmap_size = 268435456")  # 256MB in bytes https://sqlite.org/mmap.html
    cursor.execute("PRAGMA temp_store = MEMORY")    # Temporary tables in RAM
    cursor.close()


def init_huey(app):
    """
    Initialize Huey task queue configuration for the app.
    
    The actual Huey backend is configured via environment variables:
    - HUEY_BACKEND: 'memory' (default) or 'sqlite' for production
    - HUEY_SQLITE_PATH: Path to huey.db (defaults to instance/huey.db)
    - HUEY_IMMEDIATE: 'true' (default) or 'false' to run tasks async
    
    For production with background worker:
        1. Set HUEY_BACKEND=sqlite
        2. Run the Huey consumer: python -m huey.bin.huey_consumer poi_broker.worker.huey
    """
    from pathlib import Path
    
    # Ensure instance path exists (used for SQLite backend)
    Path(app.instance_path).mkdir(parents=True, exist_ok=True)
    
    # Log the active backend
    backend_name = 'SQLite' if hasattr(huey, 'filename') else 'Memory'
    if hasattr(huey, 'immediate'):
        immediate_str = ' (immediate/synchronous)' if huey.immediate else ' (async, requires worker)'
    else:
        immediate_str = ' (queue-based)'
    
    app.logger.info(f'Huey task queue initialized with {backend_name} backend{immediate_str}')

    # NOTE: stale-export cleanup is intentionally NOT run at startup. It is the
    # Huey consumer's periodic task (cleanup_stale_export_tasks, every 5 min)
    # that fails PENDING/RUNNING exports stuck past EXPORT_STALE_MAX_AGE_SECONDS.
    # Running it here would race across Gunicorn workers at boot and would hit
    # export_task before migrations/tests create it. The periodic task is the
    # single, race-free owner of this responsibility.

def create_app():
    app = Flask(__name__)
    
    """
    # Set up profiling middleware (only in development mode)
    profile_dir = "profiler_output"
    os.makedirs(profile_dir, exist_ok=True)

    # Wrap the app with ProfilerMiddleware
    app.wsgi_app = ProfilerMiddleware(
        app.wsgi_app,
        profile_dir=profile_dir,  # Save .prof files here
        restrictions=[30],        # Show top 30 functions in console
        sort_by=("cumulative",)   # Sort by cumulative time
    )
    """
    logging.basicConfig(handlers=[logging.FileHandler(filename="app.log", 
                                                 encoding='utf-8', mode='a+')],
                    format="%(asctime)s %(name)s:%(levelname)s:%(message)s", 
                    level=logging.INFO)
    # Reduce werkzeug noise
    logging.getLogger("werkzeug").setLevel(logging.ERROR)  # or logging.WARNING
    # Silence Huey's verbose debug logs (scheduler, consumer, etc.)
    logging.getLogger("huey").setLevel(logging.INFO)
    logging.getLogger("huey.consumer").setLevel(logging.INFO)
    logging.getLogger("huey.consumer.Scheduler").setLevel(logging.INFO)

    base_dir = Path(__file__).resolve().parent
    config, db_path, login_db_path = build_app_config(base_dir)
    app.config.update(config)
    app.logger.info('Configured alerts database at %s', db_path)
    app.logger.info('Configured users database at %s', login_db_path)

    if app.debug is True:
        app.jinja_env.auto_reload = True
    else:
        # Trust one hop of X-Forwarded-For/Proto/Host. x_host=1 is safe only when
        # nginx sets X-Forwarded-Host to $server_name (not the client Host).
        # Emailed links should still use PUBLIC_BASE_URL so a direct Gunicorn
        # request cannot poison verification/reset URLs.
        from werkzeug.middleware.proxy_fix import ProxyFix
        app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1)
        """
        NOTE: Make sure the NGINX site config passes the correct headers:
        location / {
            proxy_pass http://127.0.0.1:8000;
            proxy_set_header Host $host;
            proxy_set_header X-Real-IP $remote_addr;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header X-Forwarded-Proto $scheme;
            proxy_set_header X-Forwarded-Host $server_name;
        }
        """
    
    # Initialize extensions with app
    db.init_app(app)
    
    # Configure SQLite pragmas - must be inside app context to access db.engine
    with app.app_context():
        event.listens_for(db.engine, "connect")(_configure_sqlite_pragmas)
        # Also configure the users database if it's SQLite
        if 'users' in db.engines:
            event.listens_for(db.engines['users'], "connect")(_configure_sqlite_pragmas)
        _ensure_email_verification_expiry_column(app)

    init_huey(app)

    if not app.debug and not app.config.get('TESTING'):
        storage_uri = app.config.get('RATELIMIT_STORAGE_URI', 'memory://')
        if storage_uri.startswith('memory:'):
            app.logger.warning(
                'RATELIMIT_STORAGE_URI is memory://; counters are per-process. '
                'Under multi-worker Gunicorn set a shared backend URI.'
            )
        if not app.config.get('PUBLIC_BASE_URL'):
            app.logger.warning(
                'PUBLIC_BASE_URL is unset; verification and reset emails will '
                'use the request Host / X-Forwarded-Host. Set PUBLIC_BASE_URL '
                'to the canonical public origin.'
            )
    
    csrf.init_app(app)
    limiter.init_app(app)
    login_manager.init_app(app)

    # import and register blueprints here to avoid circular imports
    # Alias avoids shadowing the local Flask `app` variable (Pylance).
    from .app import register_blueprints as _register_blueprints
    from .auth import auth_blueprint
    from .classification import classification_blueprint
    from .observing_tool import observing_tool_blueprint

    _register_blueprints(app)
    app.register_blueprint(auth_blueprint)
    app.register_blueprint(observing_tool_blueprint)
    app.register_blueprint(classification_blueprint)

    # Configure Flask-Login after blueprint registration
    from .models import User
    @login_manager.user_loader
    def load_user(user_id):
        return db.session.get(User, int(user_id))
    
    # Set login_view AFTER blueprint registration to ensure the endpoint exists
    login_manager.login_view = 'auth.login'
    
    @login_manager.unauthorized_handler
    def unauthorized():
        message = 'Log-in or Sign-Up to use this feature.'
        if request.path.startswith('/api/'):
            return jsonify({'error': 'authentication required', 'message': message}), 401
        flash(message, 'warning')
        return redirect(url_for('auth.login'))

    @app.context_processor
    def inject_template_globals():
        return {
            'app_version': app.config.get('APP_VERSION', ''),
            'current_year': datetime.now(timezone.utc).year,
        }
    
    @app.after_request
    def add_security_headers(resp):
        # Enable only when behind HTTPS
        if request.is_secure:
            resp.headers.setdefault(
                'Strict-Transport-Security',
                'max-age=31536000; includeSubDomains; preload'
            )
        resp.headers.setdefault('X-Content-Type-Options', 'nosniff')
        resp.headers.setdefault('X-Frame-Options', 'DENY')
        resp.headers.setdefault('Referrer-Policy', 'strict-origin-when-cross-origin')
        resp.headers.setdefault('Permissions-Policy', 'geolocation=(), microphone=(), camera=()')
        
        resp.headers.setdefault(
            'Content-Security-Policy',
            "default-src 'self'; "
            "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://code.jquery.com https://cdnjs.cloudflare.com https://cdn.jsdelivr.net; "
            "style-src 'self' 'unsafe-inline' https://code.jquery.com https://fonts.googleapis.com https://cdn.jsdelivr.net; "
            "font-src 'self' https://fonts.gstatic.com data:; "
            "connect-src 'self' https://cdnjs.cloudflare.com https://cdn.jsdelivr.net https://storage.googleapis.com "
            "https://alaskybis.cds.unistra.fr https://aladin.cds.unistra.fr https://alasky.unistra.fr https://alasky.cds.unistra.fr https://simbad.cds.unistra.fr "
            "https://alaskybis.unistra.fr https://casda.csiro.au https://irsa.ipac.caltech.edu https://healpix.ias.u-psud.fr https://skies.esac.esa.int data:; "
            "img-src 'self' data: https:; object-src 'none'; base-uri 'self'; frame-ancestors 'none'"
        )
        return resp
    
    # Global error handler for CSRF errors raised by Flask-WTF
    @app.errorhandler(CSRFError)
    def handle_csrf_error(error):
        if request.path.startswith('/api/'):
            return jsonify({'error': 'csrf validation failed', 'message': str(error)}), 400
        flash('Your form session expired or is invalid. Please reload this page and submit again. If this is a reset link, request a new one.', 'danger')
        return redirect(request.path)

    from .cli import register_cli
    register_cli(app)

    return app
