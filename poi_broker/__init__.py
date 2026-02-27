import os
from pathlib import Path
from datetime import datetime, timedelta, UTC
import logging
from functools import lru_cache
from flask import (Flask, app, render_template, abort, jsonify, request, Response,
                   redirect, url_for, make_response, Blueprint, flash)
from flask_login import LoginManager, login_required, current_user
from flask_wtf.csrf import CSRFProtect, CSRFError
from werkzeug.middleware.proxy_fix import ProxyFix
#import jinja2
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import inspect
#import orm
from astropy.time import Time

# Initialize SQLAlchemy instance (outside create_app for import access)
db = SQLAlchemy()
login_manager = LoginManager()
csrf = CSRFProtect()

@lru_cache(maxsize=8192)
def _format_mjd_cached(mjd_value: float) -> str:
    jdate = mjd_value + 2400000.5
    dt = Time(jdate, format='jd', scale='utc').to_datetime(timezone=UTC) 
    #t = Time(jdate, format='jd')
    #dt = datetime.fromisoformat(t.isot) # ISO 8601 string in UTC
    return  dt.strftime('%Y-%m-%d %H:%M:%S')


def create_app():
    app = Flask(__name__)
    #app.config.from_pyfile(config_filename)
    app.jinja_env.auto_reload = True
    
    # Enable ProxyFix to trust headers from NGINX reverse proxy
    # This ensures url_for(..., _external=True) uses the correct X-Forwarded-* headers
    app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1)
    """
    Make sure your NGINX site config passes the correct headers:
    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_set_header X-Forwarded-Host $server_name;
    }
    """

    logging.basicConfig(handlers=[logging.FileHandler(filename="app.log", 
                                                 encoding='utf-8', mode='a+')],
                    format="%(asctime)s %(name)s:%(levelname)s:%(message)s", 
                    level=logging.ERROR)

    # Main DB (alerts)
    base_dir = Path(__file__).resolve().parent  # poi_broker directory
    # Search upward from the parent package folder
    db_path = base_dir.parent.parent / '_broker_db/ztf_alerts_stream.db'
    #print(db_path)
    db_uri = 'sqlite:///{}'.format(db_path)

    # Separate DB for logins
    login_db_path = base_dir.parent.parent / '_broker_db/users.db'
    #print(login_db_path)
    login_db_uri = 'sqlite:///{}'.format(login_db_path)

    secret_key = os.environ.get('SECRET_KEY')
    if not secret_key:
        raise RuntimeError("SECRET_KEY environment variable must be set (generate with: python -c 'import secrets; print(secrets.token_hex(32))')")

    # Load config from environment variables; use safe defaults for development
    app.config.update(
        DEBUG=os.environ.get('FLASK_DEBUG', 'False') == 'True',
        TESTING=os.environ.get('FLASK_TESTING', 'False') == 'True',
        TEMPLATES_AUTO_RELOAD=True,
        SQLALCHEMY_DATABASE_URI=db_uri,
        SQLALCHEMY_BINDS={
            'users': login_db_uri         # second DB for logins
        },
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
        SECRET_KEY=secret_key, #Used by Flask to sign the session cookie data so users can't tamper with it.
        PERMANENT_SESSION_LIFETIME = timedelta(hours=24),
        SESSION_COOKIE_SECURE = True,  # HTTPS only
        SESSION_COOKIE_HTTPONLY = True,  # Prevent XSS
        SESSION_COOKIE_SAMESITE = 'Lax'  # CSRF protection
    )
    
    # Initialize extensions with app
    db.init_app(app)
    csrf.init_app(app)
    login_manager.init_app(app)

    # import and register blueprints here to avoid circular imports
    from .app import main_blueprint
    from .observing_tool import observing_tool_blueprint
    from .classification import classification_blueprint
    from .auth import auth_blueprint

    app.register_blueprint(main_blueprint)
    app.register_blueprint(auth_blueprint)
    app.register_blueprint(observing_tool_blueprint)
    app.register_blueprint(classification_blueprint)

    # Configure Flask-Login after blueprint registration
    from .models import User
    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))
    
    # Set login_view AFTER blueprint registration to ensure the endpoint exists
    login_manager.login_view = 'auth.login'

    # Register Jinja filters
    @app.template_filter('astro_filter')
    def astro_filter(str):
        if (str == "g"):
            return "g"
        elif (str == "R"): #TODO?
            return "R"
        elif (str == "i"):
            return "i"
        else:
            return ""

    @app.template_filter('mag_filter')
    def mag_filter(num):
        if num: 
            return round(num,3)
        # else:
        #     return ''

    @app.template_filter('format_mjd_readable')
    def format_mjd_readable(value):
        if value is None:
            return ''
        
        try:
            mjd_value = float(value)
            return _format_mjd_cached(mjd_value, True) # Use Astropy for accurate conversion?
        except (TypeError, ValueError, OverflowError):
            return ''
    
    @app.errorhandler(CSRFError)
    def handle_csrf_error(error):
        flash('Your form session expired or is invalid. Please reload this page and submit again. If this is a reset link, request a new one.', 'danger')
        if request.referrer:
            return redirect(request.referrer)
        return redirect(url_for('auth.login'))
    
    return app