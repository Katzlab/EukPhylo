#!/bin/bash

# Default mode
MODE="transcriptomes"

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --mode) MODE="$2"; shift ;;  # Read the next argument as the mode
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Run the appropriate script based on the mode
if [[ "$MODE" == "transcriptomes" ]]; then
    exec /EukPhylo/PTL1/Transcriptomes/Scripts/wrapper_submit.sh
elif [[ "$MODE" == "genomes" ]]; then
    exec /EukPhylo/PTL1/Genomes/Scripts/wrapper_submit.sh
else
    echo "Invalid mode: $MODE. Use --mode transcriptomes or --mode genomes."
    exit 1
fi
