#!/bin/bash

CACHE_DIR="$HOME/.cache/openproblems/figshare/SpatialSimBench_dataset/26054188"
DATASET_ID="spatialsimbench_mobnew"
OUTPUT_DIR="resources_test/$DATASET_ID"

if [[ ! -f "$CACHE_DIR/PDAC.rds" ]]; then
  echo "Run 'scripts/create_datasets/cache.sh' first!"
  exit 1
fi

cp "$CACHE_DIR/MOBNEW.rds" "$OUTPUT_DIR/MOBNEW.rds"
cp "$CACHE_DIR/MOBNEW_sc.rds" "$OUTPUT_DIR/MOBNEW_sc.rds"

nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  -profile docker \
  --input_sc "$OUTPUT_DIR/MOBNEW_sc.rds" \
  --input_sp "$OUTPUT_DIR/MOBNEW.rds" \
  --id "$DATASET_ID" \
  --dataset_id "$DATASET_ID" \
  --dataset_name "MOBNEW" \
  --dataset_description_spatial "..." \
  --dataset_description_singlecell "..." \
  --dataset_summary_spatial "..." \
  --dataset_summary_singlecell "..." \
  --dataset_reference_spatial "..." \
  --dataset_reference_singlecell "..." \
  --dataset_url_spatial "..." \
  --dataset_url_singlecell "..." \
  --dataset_organism "mus_musculus" \
  --dataset_assay_spatial Visium \
  --dataset_assay_singlecell Chromium \
  --output_sc '$id/dataset_sc.h5ad' \
  --output_sp '$id/dataset_sp.h5ad' \
  --output_state '$id/state.yaml' \
  --publish_dir resources_test

# run a method
viash run src/methods/scdesign3/config.vsh.yaml -- \
  --input "$OUTPUT_DIR/dataset_sp.h5ad" \
  --output "$OUTPUT_DIR/simulated_dataset.h5ad"

# process the simulated spatial dataset
viash run src/process_datasets/generate_sim_spatialcluster/config.vsh.yaml -- \
  --input_sp "$OUTPUT_DIR/dataset_sp.h5ad" \
  --input_sp_sim "$OUTPUT_DIR/simulated_dataset.h5ad" \
  --output_sp "$OUTPUT_DIR/simulated_dataset_processed.h5ad"

# run a metric
viash run src/metrics/ks_statistic_spatial/config.vsh.yaml -- \
  --input_spatial_dataset "$OUTPUT_DIR/dataset_sp.h5ad" \
  --input_singlecell_dataset "$OUTPUT_DIR/dataset_sc.h5ad" \
  --input_simulated_dataset "$OUTPUT_DIR/simulated_dataset_processed.h5ad" \
  --output "$OUTPUT_DIR/score.h5ad"

# sync the resources
aws s3 sync --profile op \
  resources_test \
  s3://openproblems-data/resources_test/task_spatial_simulators \
  --delete --dryrun
