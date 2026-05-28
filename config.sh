#!/bin/bash
# Central configuration file for Thenus australiensis sex marker workflow

# Directories
export BASE_DIR="$PWD"
export DATA_DIR="$BASE_DIR/data"
export RESULT_DIR="$BASE_DIR/result"
export LOG_DIR="$BASE_DIR/logs"

# Ensure directories exist
mkdir -p "$DATA_DIR" "$RESULT_DIR" "$LOG_DIR"

# Compute resources
export THREADS=16
export MAX_MEMORY="80G"
