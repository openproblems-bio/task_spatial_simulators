#!/bin/bash

set -e

docker run --rm -it \
  -v $(pwd)/resources/task_spatial_simulators:/resources/task_spatial_simulators \
  amazon/aws-cli s3 sync \
  s3://openproblems-data/resources/task_spatial_simulators \
  /resources/task_spatial_simulators \
  --no-sign-request --delete
