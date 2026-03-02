#from app import app
from poi_broker import create_app

# create the WSGI app instance via AppFactory specified in __init__.py
app = create_app()

if __name__ == "__main__":
    app.run()
