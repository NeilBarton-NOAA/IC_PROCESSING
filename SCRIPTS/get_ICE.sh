#!/bin/bash
set -xu
dtg=${1}
SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh
source ${SCRIPT_DIR}/defaults.sh
dir=${dir_restart_ice}
mkdir -p ${dir} && cd ${dir}

echo "DOWNLOADING CICE data to ${dir}"
file_in="${aws_path}/${dtg}/gdas.${dtg:0:8}/${dtg:8:10}/model/ice/restart/${DTG_TEXT_SRC}.cice_model.res.nc"
file_out=${DTG_TEXT_DES}.cice_model.res.nc

WGET_AWS ${file_in} ${file_out} 
FIND_EMPTY_FILES ${PWD}

exit 0
