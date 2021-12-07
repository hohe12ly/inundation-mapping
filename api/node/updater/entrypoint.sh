#!/bin/sh

umask 002
cd /opt/updater/
echo "==============================="
echo "Starting Update Loop"
date
python3 ./updater.py
