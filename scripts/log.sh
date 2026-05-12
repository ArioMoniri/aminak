#!/bin/bash
# Helper: log a line with timestamp to pipeline.log
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$PROJECT_DIR/pipeline.log"
}
export -f log
