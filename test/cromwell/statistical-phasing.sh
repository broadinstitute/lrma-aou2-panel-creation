#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

WORKING_DIR=/home/runner/work
WDL_DIR=$WORKING_DIR/wdl
CROMWELL_TEST_DIR=$WORKING_DIR/test/cromwell

set -e
echo "Test..."
