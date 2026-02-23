from pathlib import Path
from datetime import datetime, timedelta
from flask import (Flask, app, render_template, abort, jsonify, request, Response,
                   redirect, url_for, logging, make_response, Blueprint)
from flask_login import LoginManager, login_required, current_user
import jinja2
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import inspect
#import orm
from astropy.time import Time

# Initialize SQLAlchemy instance (outside create_app for import access)
db = SQLAlchemy()
login_manager = LoginManager()

def create_app():
    app = Flask(__name__)
    #app.config.from_pyfile(config_filename)
    app.jinja_env.auto_reload = True

    # Main DB (alerts)
    base_dir = Path(__file__).resolve().parent  # poi_broker directory
    # Search upward from the parent package folder
    db_path = base_dir.parent / '_broker_db/ztf_alerts_stream.db'
    #print(db_path)
    db_uri = 'sqlite:///{}'.format(db_path)

    # Separate DB for logins
    login_db_path = base_dir.parent / '_broker_db/users.db'
    #print(login_db_path)
    login_db_uri = 'sqlite:///{}'.format(login_db_path)


    #TODO: Move to .env o.Ae.
    app.config.update(
        DEBUG=True,
        TESTING=True,
        TEMPLATES_AUTO_RELOAD=True,
        SQLALCHEMY_DATABASE_URI=db_uri,
        SQLALCHEMY_BINDS={
            'users': login_db_uri         # second DB for logins
        },
        SQLALCHEMY_TRACK_MODIFICATIONS=False,
        SECRET_KEY = 'your-secret-key-change-in-production',
        PERMANENT_SESSION_LIFETIME = timedelta(hours=24),
        SESSION_COOKIE_SECURE = True,  # HTTPS only
        SESSION_COOKIE_HTTPONLY = True,  # Prevent XSS
        SESSION_COOKIE_SAMESITE = 'Lax'  # CSRF protection
    )
    
    # Initialize extensions with app
    db.init_app(app)
    
    # Configure Flask-Login
    #login_manager = LoginManager()
    login_manager.login_view = 'auth.login'
    login_manager.init_app(app)
    
    # User loader function for Flask-Login
    from .models import User
    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))

    #Jinja Filter
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
        #return (value / 86400000) + 40587
        jdate = value+2400000.5
        t = Time(jdate, format='jd')
        #return t.isot #in UTC
        dt = datetime.fromisoformat(t.isot) # ISO 8601 string in UTC
        return  dt.strftime('%Y-%m-%d %H:%M:%S')

    # import and register blueprints here to avoid circular imports
    from .app import main_blueprint
    from .observing_tool import observing_tool_blueprint
    from .classification import classification_blueprint
    from .auth import auth_blueprint

    app.register_blueprint(main_blueprint)
    app.register_blueprint(auth_blueprint)
    app.register_blueprint(observing_tool_blueprint)
    app.register_blueprint(classification_blueprint)

        
    return app