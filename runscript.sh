#!/bin/bash
# Host-provided arguments
HOST_INPUT="$1"
HOST_OUTPUT="$2"


# Bound container paths
IN_MOUNT="/mnt/in"
OUT_MOUNT="/mnt/out"

# File name extraction
INPUT_FILENAME=$(basename "$HOST_INPUT")
CONTAINER_INPUT="${IN_MOUNT}/${INPUT_FILENAME}"
BASENAME="${INPUT_FILENAME%.*}"

cd /opt/icefinder2
exec /opt/micromamba/bin/micromamba run -p /opt/micromamba/envs/icefinder2_env python ICEfinder2.py -i $CONTAINER_INPUT -t Single &
wait

# Expected output directory
ICEFINDER_RESULT="/opt/icefinder2/result/$BASENAME"

# Target location for the copied result
COPY_DEST="${OUT_MOUNT}/${BASENAME}"

# Ensure the output directory exists
mkdir -p "$OUT_MOUNT"

# Copy result if it exists
if [ -d "$ICEFINDER_RESULT" ]; then
    cp -r "$ICEFINDER_RESULT" "$COPY_DEST"
else
    echo "ICEfinder failed or result not found at $ICEFINDER_RESULT" >&2
    exit 1
fi
