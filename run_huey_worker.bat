@echo off
REM Huey consumer startup script for production (Windows)
REM
REM This script starts the Huey worker process that processes background tasks.
REM
REM Usage:
REM     run_huey_worker.bat
REM
REM Requirements:
REM     - Flask app environment must be configured (SECRET_KEY, DB paths, etc.)
REM     - HUEY_BACKEND must be set to 'sqlite' in the environment
REM
REM Example full production setup:
REM     set SECRET_KEY=your-secret-key
REM     set HUEY_BACKEND=sqlite
REM     set HUEY_SQLITE_PATH=C:\var\lib\poi_broker\huey.db
REM     set ALERTS_DB_PATH=C:\var\lib\poi_broker\ztf_alerts_stream.db
REM     set USERS_DB_PATH=C:\var\lib\poi_broker\users.db
REM     
REM     # Run the Flask app in one terminal
REM     python -m gunicorn -w 4 -b 0.0.0.0:8000 wsgi:app
REM     
REM     # Run Huey worker in another terminal
REM     run_huey_worker.bat

setlocal enabledelayedexpansion

REM Check if HUEY_BACKEND is set to sqlite
if not "%HUEY_BACKEND%"=="sqlite" (
    echo Warning: HUEY_BACKEND is not set to 'sqlite'
    echo Set HUEY_BACKEND=sqlite for async task processing
    echo Without this, tasks will run synchronously (development mode)
    echo.
)

REM Determine number of worker threads
if "%HUEY_THREADS%"=="" (
    set THREADS=4
) else (
    set THREADS=%HUEY_THREADS%
)

REM Determine logfile location
if "%HUEY_LOGFILE%"=="" (
    set LOGFILE=huey.log
) else (
    set LOGFILE=%HUEY_LOGFILE%
)

echo Starting Huey consumer...
echo   Backend: %HUEY_BACKEND:~0,6%
if not "%HUEY_SQLITE_PATH%"=="" (
    echo   Database: %HUEY_SQLITE_PATH%
)
echo   Worker threads: %THREADS%
echo   Logfile: %LOGFILE%
echo.

REM Run the Huey consumer
REM See https://huey.readthedocs.io/en/latest/cli.html for CLI options
python -m huey_consumer ^
    poi_broker.extensions.huey ^
    --workers=%THREADS% ^
    --worker-type=thread ^
    --logfile=%LOGFILE% ^
    --verbose

endlocal
