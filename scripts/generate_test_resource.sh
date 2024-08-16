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
  --publish_dir resources_test/datasets

viash run src/methods/scdesign3/config.vsh.yaml -- \
  --input resources_test/datasets/MOBNEW/dataset_sp.h5ad \
  --output resources_test/datasets/MOBNEW/simulated_dataset.h5ad

viash run src/metrics/ks_statistic/config.vsh.yaml -- \
  --input_spatial_dataset resources_test/datasets/MOBNEW/dataset_sp.h5ad \
  --input_singlecell_dataset resources_test/datasets/MOBNEW/dataset_sc.h5ad \
  --input_simulated_dataset resources_test/datasets/MOBNEW/simulated_dataset.h5ad \
  --output resources_test/datasets/MOBNEW/score.h5ad
