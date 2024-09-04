#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# NOTE: depending on the the datasets and components, you may need to launch this workflow
# on a different compute platform (e.g. a HPC, AWS Cloud, Azure Cloud, Google Cloud).
# please refer to the nextflow information for more details:
# https://www.nextflow.io/docs/latest/

set -e

# generate a unique id
RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="resources/results/${RUN_ID}"

# write the parameters to file
cat > /tmp/params.yaml << HERE
input_states: resources/datasets/**/state.yaml
rename_keys: 'input_singlecell_dataset:output_sc;input_spatial_dataset:output_sp'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# run the benchmark
nextflow run openproblems-bio/task_spatial_simulators \
  --revision build/main \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c common/nextflow_helpers/labels.config \
  -params-file /tmp/params.yaml
