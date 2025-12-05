#!/bin/bash -l

set -euxo pipefail

# cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

# set repo dir if we are not running locally and have already set it, e.g.:
# sudo REPO_DIR=/local/root-of-this-repo CROMWELL_JAR=/local/cromwell.jar ./this-script.sh
REPO_DIR=${REPO_DIR:=/Users/suhang/Analysis/lrma-aou2-panel-creation/}

# insert repo dir into resource files (this will create *.mod.* files, which may need to be cleaned up locally)
sed -e "s|__REPO_DIR__|$REPO_DIR|g" $REPO_DIR/test/resources/statistical-phasing/statistical-phasing.json > $REPO_DIR/test/resources/statistical-phasing/statistical-phasing.mod.json
sed -e "s|__REPO_DIR__|$REPO_DIR|g" $REPO_DIR/test/resources/statistical-phasing/genetic_map_b38.tsv > $REPO_DIR/test/resources/statistical-phasing/genetic_map_b38.mod.tsv

java -jar $CROMWELL_JAR run $REPO_DIR/wdl/methods/phasing/StatisticalPhasing_phase2.wdl -i $REPO_DIR/test/resources/statistical-phasing/statistical-phasing.mod.json
