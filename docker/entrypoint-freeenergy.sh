#!/bin/bash
set -e

MODEL_PATH="${NEQUIP_MODEL_PATH:-/opt/nequip/NequIP-OAM-XL.nequip.pt2}"

if [ ! -f "$MODEL_PATH" ]; then
    echo "Compiling NequIP-OAM-XL model (first run only)..."
    mkdir -p "$(dirname "$MODEL_PATH")"
    nequip-compile \
        nequip.net:mir-group/NequIP-OAM-XL:0.1 \
        "$MODEL_PATH" \
        --device cuda \
        --mode aotinductor \
        --target ase
    echo "Model compiled to $MODEL_PATH"
fi

exec "$@"
