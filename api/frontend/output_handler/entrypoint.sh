#!/bin/sh

umask 002
cd /opt/output_handler/
echo "==============================="
echo "Starting Output Handler"
date
python ./output_handler.py