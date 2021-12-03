#!/bin/sh

cd /opt/gui/
echo "==============================="
echo "Starting Gunicorn"
date
exec gunicorn --bind 0.0.0.0:5000 --reload wsgi:app