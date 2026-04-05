from pathlib import Path
from datetime import datetime, timezone
import logging
from flask import (Flask, render_template, abort, jsonify, request, Response,
                   redirect, url_for, make_response, Blueprint, flash)
from flask_login import LoginManager, login_required, current_user
from flask_wtf.csrf import CSRFProtect, CSRFError
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from sqlalchemy import event

#from werkzeug.middleware.profiler import ProfilerMiddleware
#import jinja2
from flask_sqlalchemy import SQLAlchemy
from astropy.time import Time
from .settings import build_app_config

# Initialize SQLAlchemy instance (outside create_app for import access)
db = SQLAlchemy()
login_manager = LoginManager()
csrf = CSRFProtect()
limiter = Limiter(key_func=get_remote_address, default_limits=[])

def _configure_sqlite_pragmas(dbapi_conn, connection_record):
    """Configure SQLite PRAGMAs for improved concurrent read performance."""
    cursor = dbapi_conn.cursor()
    cursor.execute("PRAGMA journal_mode = WAL") # Enable Write-Ahead Logging for better concurrency https://sqlite.org/wal.html
    cursor.execute("PRAGMA synchronous = NORMAL") # Balance durability and performance; WAL mode is safe with NORMAL https://sqlite.org/wal.html#durability
    cursor.execute("PRAGMA cache_size = -64000")  # 64MB in KB https://sqlite.org/pragma.html#pragma_cache_size
    cursor.execute("PRAGMA mmap_size = 268435456")  # 256MB in bytes https://sqlite.org/mmap.html
    cursor.execute("PRAGMA temp_store = MEMORY")    # Temporary tables in RAM
    cursor.close()

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

    base_dir = Path(__file__).resolve().parent
    config, db_path, login_db_path = build_app_config(base_dir)
    app.config.update(config)
    app.logger.info('Configured alerts database at %s', db_path)
    app.logger.info('Configured users database at %s', login_db_path)

    if app.debug is True:
        app.jinja_env.auto_reload = True
    else:
        # Enable ProxyFix to trust headers from NGINX reverse proxy in production
        from werkzeug.middleware.proxy_fix import ProxyFix
        # This ensures url_for(..., _external=True) uses the correct X-Forwarded-* headers
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
    
    csrf.init_app(app)
    limiter.init_app(app)
    login_manager.init_app(app)

    # import and register blueprints here to avoid circular imports
    from .app import register_blueprints
    from .observing_tool import observing_tool_blueprint
    from .classification import classification_blueprint
    from .auth import auth_blueprint

    register_blueprints(app)
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
    
    # Global error handler for CSRF errors raised by Flask-WTF
    @app.errorhandler(CSRFError)
    def handle_csrf_error(error):
        if request.path.startswith('/api/'):
            return jsonify({'error': 'csrf validation failed', 'message': str(error)}), 400
        flash('Your form session expired or is invalid. Please reload this page and submit again. If this is a reset link, request a new one.', 'danger')
        if request.referrer:
            return redirect(request.referrer)
        return redirect(url_for('auth.login'))
    
    return app