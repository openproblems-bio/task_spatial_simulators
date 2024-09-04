#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# remove this when you have implemented the script
echo "TODO: once the 'run_benchmark' workflow has been implemented, update this script to use it."
echo "  Step 1: replace 'task_spatial_simulators' with the name of the task in the following command."
echo "  Step 2: replace the rename keys parameters to fit your run_benchmark inputs"
echo "  Step 3: replace the settings parameter to fit your run_benchmark outputs"
echo "  Step 4: remove this message"
exit 1

set -e

# write the parameters to file
cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources_test/task_spatial_simulators/**/state.yaml
rename_keys: 'input_singlecell_dataset:output_sc;input_spatial_dataset:output_sp'
output_state: "state.yaml"
publish_dir: s3://openproblems-nextflow/temp/task_spatial_simulators/
HERE

nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c common/nextflow_helpers/labels.config \
  -params-file /tmp/params.yaml
