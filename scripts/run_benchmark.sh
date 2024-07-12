#!/bin/bash

OUTPUT_DIR=output/tmp
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

# nextflow run . \
#   -main-script target/nextflow/workflows/run_benchmark/main.nf \
#   -profile docker \
#   -resume \
#   -entry auto \
#   -c common/nextflow_helpers/labels_ci.config \
#   --input_states "resources/datasets/**/state.yaml" \
#   --rename_keys 'input_singlecell_dataset:output_sc;input_spatial_dataset:output_sp' \
#   --settings '{"output_scores": "scores.yaml", "output_dataset_info": "dataset_info.yaml", "output_method_configs": "method_configs.yaml", "output_metric_configs": "metric_configs.yaml", "output_task_info": "task_info.yaml"}' \
#   --publish_dir $OUTPUT_DIR \
#   --output_state "state.yaml"

nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -c common/nextflow_helpers/labels_ci.config \
  --id test \
  --input_spatial_dataset resources_test/datasets/MOBNEW/dataset_sp.h5ad \
  --input_singlecell_dataset resources_test/datasets/MOBNEW/dataset_sc.h5ad \
  --publish_dir $OUTPUT_DIR \
  --output_state "state.yaml"