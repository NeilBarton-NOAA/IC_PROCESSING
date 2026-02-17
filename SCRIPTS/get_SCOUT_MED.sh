#!/bin/bash
set -u
dtg=${1}

SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh
source ${SCRIPT_DIR}/defaults.sh
dir=${dir_restart_med}
mkdir -p ${dir} && cd ${dir}
echo "DOWNLOADING MEDIATOR data to ${dir}"

f=${DTG_TEXT_SRC}.ufs.cpld.cpl.r.nc
file_in=${aws_path}/${dtg}/gdas.${dtg:0:8}/${dtg:8:10}/model/med/restart/${f}
file_out=${DTG_TEXT_DES}.ufs.cpld.cpl.r.nc
WGET_AWS ${file_in} ${file_out}

FIND_EMPTY_FILES ${PWD}

exit 0
