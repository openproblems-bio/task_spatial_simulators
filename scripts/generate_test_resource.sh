#!/bin/bash

nextflow run . \
  -main-script target/nextflow/workflows/process_datasets/main.nf \
  --input_sc resources_test/datasets_raw/MOBNEW/dataset_sc.rds \
  --input_sp resources_test/datasets_raw/MOBNEW/dataset_sp.rds \
  --output_sc resources_test/datasets/MOBNEW/dataset_sc.h5ad \
  --output_sp resources_test/datasets/MOBNEW/temp_dataset_sp_part1.h5ad \
  --id MOBNEW \
  --dataset_id MOBNEW \
  --dataset_name "MOBNEW" \
  --dataset_description "MOBNEW" \
  --dataset_summary "MOBNEW" \
  --dataset_reference "..." \
  --dataset_organism "mus_musculus" \
  --output_state '$id/state.yaml' \
  --publish_dir resources_test/datasets

# run a method
viash run src/methods/scdesign3/config.vsh.yaml -- \
  --input resources_test/datasets/MOBNEW/dataset_sp.h5ad \
  --output resources_test/datasets/MOBNEW/simulated_dataset.h5ad

# process the simulated spatial dataset
viash run src/process_datasets/generate_sim_spatialcluster/config.vsh.yaml -- \
  --input_sp resources_test/datasets/MOBNEW/dataset_sp.h5ad \
  --input_sp_sim resources_test/datasets/MOBNEW/simulated_dataset.h5ad \
  --output_sp resources_test/datasets/MOBNEW/simulated_dataset_processed.h5ad

# run a metric
viash run src/metrics/ks_statistic_spatial/config.vsh.yaml -- \
  --input_spatial_dataset resources_test/datasets/MOBNEW/dataset_sp.h5ad \
  --input_singlecell_dataset resources_test/datasets/MOBNEW/dataset_sc.h5ad \
  --input_simulated_dataset resources_test/datasets/MOBNEW/simulated_dataset_processed.h5ad \
  --output resources_test/datasets/MOBNEW/score.h5ad

# sync the resources
aws s3 sync --profile op \
  resources_test \
  s3://openproblems-data/resources_test/task_spatial_simulators \
  --delete --dryrun
