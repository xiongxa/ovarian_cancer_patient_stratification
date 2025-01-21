#!/bin/bash

# Gets the directory where the script is currently located
PARP_DIR="$(cd "$(dirname "$0")" && pwd)"

# Prints the directory where the script is currently located
echo "PWD: $PARP_DIR"
export PARP_DIR
export PYTHONPATH="$PYTHONPATH:$PARP_DIR"
