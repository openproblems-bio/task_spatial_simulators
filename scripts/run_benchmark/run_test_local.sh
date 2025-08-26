#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# remove this when you have implemented the script
# echo "TODO: once the 'run_benchmark' workflow has been implemented, update this script to use it."
# echo "  Step 1: replace 'task_spatial_simulators' with the name of the task in the following command."
# echo "  Step 2: replace the rename keys parameters to fit your run_benchmark inputs"
# echo "  Step 3: replace the settings parameter to fit your run_benchmark outputs"
# echo "  Step 4: remove this message"
# exit 1

set -e

# write the parameters to file
cat > /tmp/params.yaml << 'HERE'
input_states: "resources_test/spatialsimbench_mobnew/state.yaml"
rename_keys: "input_singlecell_dataset:output_sc;input_spatial_dataset:output_sp"
settings: '{"output_scores": "scores.yaml", "output_dataset_info": "dataset_info.yaml", "output_method_configs": "method_configs.yaml", "output_metric_configs": "metric_configs.yaml", "output_task_info": "task_info.yaml"}'
output_state: "state.yaml"
publish_dir: "output"
HERE

nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file /tmp/params.yaml
