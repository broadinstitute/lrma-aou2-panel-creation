#!/bin/bash -l

set -euxo pipefail

# cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

# set repo dir if we are not running locally and have already set it, e.g.:
# sudo REPO_DIR=/local/root-of-this-repo CROMWELL_JAR=/local/cromwell.jar ./this-script.sh
# to run this correctly locally, you need java 17, cromwell jar, git-lfs and bcftools installed
# and you need to set jq path if it is not in your PATH
REPO_DIR=${REPO_DIR:=/home/runner/work/lrma-aou2-panel-creation/lrma-aou2-panel-creation}

# insert repo dir into resource files (this will create *.mod.* files, which may need to be cleaned up locally)
sed -e "s|__REPO_DIR__|$REPO_DIR|g" $REPO_DIR/test/resources/pangenie-panel-creation/pangenie-panel-creation.json > $REPO_DIR/test/resources/pangenie-panel-creation/pangenie-panel-creation.mod.json

java -jar $CROMWELL_JAR run $REPO_DIR/wdl/methods/imputation/PanGeniePanelCreation.wdl -i $REPO_DIR/test/resources/pangenie-panel-creation/pangenie-panel-creation.mod.json -m $REPO_DIR/test/resources/pangenie-panel-creation/pangenie-panel-creation.mod.output.json

RESULT=$(jq -r '.outputs."PanGeniePanelCreation.panel_vcf_gz"' $REPO_DIR/test/resources/pangenie-panel-creation/pangenie-panel-creation.mod.output.json)
EXPECTED=$REPO_DIR/test/resources/large/pangenie-panel-creation/expected/40-HPRC-1kGP.chr6-70M-80M.phased.ligated.snippeted-1M.shifted.prepare.id.split.mergehap.vcf.gz

diff <(bcftools view -h --no-version $EXPECTED | grep -v fileDate) <(bcftools view -h --no-version $RESULT | grep -v fileDate)
diff <(bcftools view -H --no-version $EXPECTED | grep -v fileDate) <(bcftools view -H --no-version $RESULT | grep -v fileDate) | head -500
