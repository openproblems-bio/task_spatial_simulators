#!/bin/bash

viash run src/process_datasets/convert/config.vsh.yaml -- \
  --input_sc resources_test/datasets_raw/MOBNEW/dataset_sc.rds \
  --input_sp resources_test/datasets_raw/MOBNEW/dataset_sp.rds \
  --output_sc resources_test/datasets/MOBNEW/dataset_sc.h5ad \
  --output_sp resources_test/datasets/MOBNEW/dataset_sp.h5ad \
  --dataset_id MOBNEW \
  --dataset_name "MOBNEW" \
  --dataset_description "MOBNEW" \
  --dataset_summary "MOBNEW" \
  --dataset_reference "..." \
  --dataset_organism "mus_musculus"
