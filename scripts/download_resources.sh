#!/bin/bash

set -e

echo ">> Downloading resources"

common/sync_resources/sync_resources \
  --input "s3://openproblems-data/resources_test/spatial_simulators/" \
  --output "resources_test" \
  --delete

# aws s3 sync --profile op resources_test s3://openproblems-data/resources_test/task_spatial_simulators/ --dryrun