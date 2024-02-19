#!/bin/bash
set -x

find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} \;
find . -type d -name "__pycache__" -exec rm -rf {} \;
find . -name ".DS_Store" -delete
