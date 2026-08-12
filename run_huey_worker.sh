#!/usr/bin/env bash
#
# Huey consumer startup script for production
# 
# This script starts the Huey worker process that processes background tasks.
# 
# Usage:
#     ./run_huey_worker.sh
# 
# Requirements:
#     - Flask app environment must be configured (SECRET_KEY, DB paths, etc.)
#     - HUEY_BACKEND must be set to 'sqlite' in the environment
#
# Example full production setup:
#     export SECRET_KEY=your-secret-key
#     export HUEY_BACKEND=sqlite
#     export HUEY_SQLITE_PATH=/var/lib/poi_broker/huey.db
#     export ALERTS_DB_PATH=/var/lib/poi_broker/ztf_alerts_stream.db
#     export USERS_DB_PATH=/var/lib/poi_broker/users.db
#     
#     # Run the Flask app in one terminal
#     gunicorn -w 4 -b 0.0.0.0:8000 wsgi:app
#     
#     # Run Huey worker in another terminal
#     ./run_huey_worker.sh
#

set -e

# Check if HUEY_BACKEND is set to sqlite
if [[ "${HUEY_BACKEND}" != "sqlite" ]]; then
    echo "Warning: HUEY_BACKEND is not set to 'sqlite'"
    echo "Set HUEY_BACKEND=sqlite for async task processing"
    echo "Without this, tasks will run synchronously (development mode)"
fi

# Determine number of worker threads
THREADS=${HUEY_THREADS:-2}

# Determine logfile location
LOGFILE=${HUEY_LOGFILE:-huey.log}

echo "Starting Huey consumer..."
echo "  Backend: ${HUEY_BACKEND:-memory}"
if [[ -n "${HUEY_SQLITE_PATH}" ]]; then
    echo "  Database: ${HUEY_SQLITE_PATH}"
fi
echo "  Worker threads: ${THREADS}"
echo "  Logfile: ${LOGFILE}"
echo ""

# Run the Huey consumer
# See https://huey.readthedocs.io/en/latest/cli.html for CLI options
# NOTE: target poi_broker.worker.huey (imports tasks so they register on the
# shared Huey instance). Pointing at poi_broker.extensions.huey would start a
# worker that never executes any task.
# The importable module is huey.bin.huey_consumer (NOT huey_consumer, which is
# only a console-script entry point in the venv bin). Using `python -m` avoids
# PATH issues.
exec python -m huey.bin.huey_consumer \
    poi_broker.worker.huey \
    --workers=${THREADS} \
    --worker-type=thread \
    --logfile=${LOGFILE} \
    --verbose
