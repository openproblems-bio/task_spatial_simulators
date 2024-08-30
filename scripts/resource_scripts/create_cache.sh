#!/bin/bash

# download and process all datasets from figshare
# https://figshare.com/articles/dataset/SpatialSimBench_dataset/26054188
# and store them in ~/.cache/openproblems/figshare/SpatialSimBench_dataset/26054188

DOWNLOAD_URL="https://figshare.com/ndownloader/articles/26054188/versions/2"
ZIP_FILE="$HOME/.cache/openproblems/figshare/SpatialSimBench_dataset/26054188.zip"
CACHE_DIR="$HOME/.cache/openproblems/figshare/SpatialSimBench_dataset/26054188"

if [[ ! -f "$ZIP_FILE" ]]; then
  wget "$DOWNLOAD_URL" -O "$ZIP_FILE"
fi

if [[ ! -f "$CACHE_DIR/PDAC.rds" ]]; then
  unzip "$ZIP_FILE" -d "$CACHE_DIR"
  mv "$CACHE_DIR/47112109_MOBNEW.rds" "$CACHE_DIR/MOBNEW.rds"
  mv "$CACHE_DIR/47113516_PDAC.rds" "$CACHE_DIR/PDAC.rds"
  rm "$CACHE_DIR/47115889_MOBNEW.rds"
  rm "$CACHE_DIR/47115904_PDAC.rds"
fi

