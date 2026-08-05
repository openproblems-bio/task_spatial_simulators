#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

DATASETS_DIR="resources_test/common"
OUTPUT_DIR="output/process_datasets_test"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

export NXF_VER=24.04.3

nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -c common/nextflow_helpers/labels_ci.config \
  --id run_test \
  --input_sc "$DATASETS_DIR/MOBNEW_sc.rds" \
  --input_sp "$DATASETS_DIR/MOBNEW.rds" \
  --dataset_id "spatialsimbench_mobnew" \
  --dataset_name "MOBNEW" \
  --dataset_description "..." \
  --dataset_summary "..." \
  --dataset_reference "..." \
  --dataset_url "..." \
  --dataset_organism "mus_musculus" \
  --dataset_assay_spatial Visium \
  --dataset_assay_singlecell Chromium \
  --output_sc "dataset_sc.h5ad" \
  --output_sp "dataset_sp.h5ad" \
  --output_state "state.yaml" \
  --publish_dir "$OUTPUT_DIR" \
  --output_state "state.yaml"
