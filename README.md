# Point of Interest Community Broker

The Point of Interest Community Broker is a transient alert software (Rubin Observatory (LSST) Alert and Community Broker) that is currently tested with the ZTF alert stream.
Incoming alerts will be processed, annotated, classified and forwarded.

For a general overview of LSST Community Broker, see https://www.lsst.org/scientists/alert-brokers


This repository contains the frontend (website) for the Point of Interest broker.
The Point of Interest broker can be found at https://poibroker.uantof.cl/


## Motivation
 
Current and especially upcoming all-sky time-domain surveys, such as LSST, will deliver a vast amount of data each night, requiring for the developent of flexible, straightforward tools for the analysis, selection and forwarding of information regarding astrophysical transients and variable objects. 
 
Our alert broker, called *Point of Interest*, is tailored towards the needs of astronomers looking for updated observations of variable stars in specific on-sky regions. Developed by a small team at Vanderbilt University, where I'm the main developer responsible for this project, this *Point of Interest*' alert broker should enable users to get updates on variable star observations from a straightforward, user-friendly web service. Data are processed in real time by big data/ machine learning algorithms and will be immediately available to the user community.


*Point of Interest* differs from other brokers in the focus on updates on variable stars, thus running a rather specific than the full analysis chain of streamed data. As a consequence, the broker is rather lightweight. *Point of Interest* users are encouraged to design their own on-sky regions they want receive updates for (such as for planned follow-up campaigns) or select from a list of on-sky regions which are particularly interesting for variable star observers, such as stellar streams, globular clusters and dwarf galaxies.

## Usage
This repository contais the web frontend, including a small database for testing purposes.

After downloading the package, within the folder `_web_frontend` install all required packages:

`pip3 install -r requirements.txt`

> [!NOTE]
> requirements2.txt shows the output of `pip freeze` and lists exact version numbers of the packages. However, this is more for reference. You may not be able to install these exact version on your specific system/environment.

Start the web app from the terminal from the **project root** (the folder that **contains** the _sub-folder_ `poi_boker`) with 

`flask --app "poi_broker:create_app()" run --debug`

or 

`flask --app "poi_broker:create_app()" run --no-debug --no-reload`

In the web browser, enter `http://127.0.0.1:5000/` to view the front-ent.

In case the website isn't displayed: do a

`cat app.log`

in your terminal window to see the correct URL


Also, inspect the browser developer console (F12) to see if there are any JavaScript errors, ie. a missmatching bokeh.min.js version.

The app can be profiled like any other python module, e.g.

`python -m cProfile -o program.prof -m flask --app "poi_broker:create_app()" run --no-reload`

but it's more useful to use the werkzeug ProfilerMiddleware to profile routes. Activate it in `__init__.py` if needed.
