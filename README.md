# Point of Interest Community Broker

The Point of Interest Community Broker is a transient alert software (Rubin Observatory (LSST) Alert and Community Broker) that is currently tested with the ZTF alert stream.
Incoming alerts will be processed, annotated, classified and forwarded.

For a general overview of LSST Community Broker, see https://www.lsst.org/scientists/alert-brokers


This repository contains the frontend (website) for the Point of Interest broker.
The Point of Interest broker can be found at https://poibroker.uantof.cl/


## Motivation
 
Current and especially upcoming all-sky time-domain surveys, such as LSST, will deliver a vast amount of data each night, requiring for the development of flexible, straightforward tools for the analysis, selection and forwarding of information regarding astrophysical transients and variable objects. 
 
Our alert broker, called *Point of Interest*, is tailored towards the needs of astronomers looking for updated observations of variable stars in specific on-sky regions. Developed by a small team at Vanderbilt University, where I'm the main developer responsible for this project, this *Point of Interest*' alert broker should enable users to get updates on variable star observations from a straightforward, user-friendly web service. Data are processed in real time by big data/ machine learning algorithms and will be immediately available to the user community.


*Point of Interest* differs from other brokers in the focus on updates on variable stars, thus running a rather specific than the full analysis chain of streamed data. As a consequence, the broker is rather lightweight. *Point of Interest* users are encouraged to design their own on-sky regions they want receive updates for (such as for planned follow-up campaigns) or select from a list of on-sky regions which are particularly interesting for variable star observers, such as stellar streams, globular clusters and dwarf galaxies.

## Installation
This repository contains the web frontend, including a small database for testing purposes.

After downloading the package, install all required packages for development:

`pip3 install -r requirements-dev.txt`

Use the regular `requirements.txt` for deployments.

> [!NOTE]
> requirements2.txt shows the output of `pip freeze` and lists exact version numbers of the packages. However, this is mostly for reference. You may not be able to install these exact version on your specific system/environment.

Create a `.env` file in the project root based on the `.env.example` file and adjust the settings as needed.

Database locations are configured through environment variables:

- `ALERTS_DB_PATH`: path to the main alerts SQLite database.
- `USERS_DB_PATH`: path to the users/auth SQLite database.

Both values may be absolute paths or paths relative to the project root. By default the app assumes both files live in a sibling directory named `_broker_db` outside the repository, for example:

```
../_broker_db/ztf_alerts_stream.db
../_broker_db/users.db
```

> [!NOTE]
> If you don't have the database files, you can create them using the alerts and user database schema files found in the `/tools` folder. Sample alerts data can be found in the same directory.

```
sqlite3 alerts.db < alertsdb_schema.sql
sqlite3 users.db < usersdb_schema.sql 
sqlite3 alerts.db < alerts_dump.sql 
```
The main alerts database is intentionally excluded because it is too large to keep inside the repository.


## Usage

Start the web app from the terminal from the **project root** (the folder that **contains** the _sub-folder_ `poi_boker`) with 

`flask --app "poi_broker:create_app()" run --debug`

or 

`flask --app "poi_broker:create_app()" run --no-debug --no-reload`

In the web browser, enter `http://127.0.0.1:5000/` to view the front-ent.

In case the website isn't displayed: do a

`cat app.log`

in your terminal window to see the correct URL.


Also, inspect the browser developer console (F12) to see if there are any errors, e.g. missing includes or JavaScript exceptions.

### Testing, Debugging and Profiling

For now, only a basic `pytest` smoke test for the main route exists. Run the tests before commiting changes.

```
(poi_brokerenv) λ pytest -q
```

It is highly recommended to debug the app in a capable IDE like VS Code to leverage built-in debugging capabilities.

To set breakpoints in the CLI use pdb
```
import pdb
from pprint import pprint
...
pdb.set_trace()
```

Use the logger class to log output to the app.log file instead of `print`'ing to the console
```
app.logger.info('Info')
app.logger.warning('Warn')
logging.error('Exception occurred', exc_info=True) #or: logging.exception()
```

The app can be profiled like any other python module, e.g.

```
python -m cProfile -o program.prof -m flask --app "poi_broker:create_app()" run --no-reload
```

However, it's more useful to use werkzeug's `ProfilerMiddleware` to profile routes. Activate it in `poi_broker/__init__.py` if needed.
