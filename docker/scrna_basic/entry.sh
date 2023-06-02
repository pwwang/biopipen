#!/usr/bin/env bash

# Check if we have any arguments
if [ $# -eq 0 ]; then
    echo "There are two ways to run the pipeline:"
    echo "1. Run the pipeline from the command line"
    echo "   docker/singularity run <options> pipen run scrna_basic ScrnaBasic @<configfile> [options]"
    echo "2. Run the pipeline using pipen-board:"
    echo "   docker/singularity run <options> pipen board biopipen.ns.scrna_basic:ScrnaBasic -a /workdir/board.toml [options]"
    exit 1
else
    # Run the command passed from the command line
    exec "$@"
fi
