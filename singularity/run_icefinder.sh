#!/bin/bash
# offer place for binds
mkdir -p /tmp/res
mkdir -p /tmp/tmp

# get apptainer ready
ml apptainer
IMAGE="$1"  # path to icefinder.sif image

# get inputs
INPUT="$2"
PARENT_DIR=$(dirname "$INPUT")
OUTPUT="$3"

# perform actual run
apptainer run \
    --bind "$PARENT_DIR":/mnt/in \
    --bind "$OUTPUT":/mnt/out \
    --bind /tmp/tmp:/opt/icefinder2/tmp \
    --bind /tmp/res:/opt/icefinder2/result \
    "$IMAGE" "$INPUT" "$OUTPUT" &
wait

# clean up
rm -rf /tmp/tmp
rm -rf /tmp/res
