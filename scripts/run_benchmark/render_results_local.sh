#!/bin/bash

# fail on error
set -e

# ensure we're in the root of the repo
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

publish_dir="resources/results/...some_run_id..."

common/scripts/render_results_report local "$publish_dir" --output "$publish_dir/report/"
