#!/bin/bash
set -u
export dtg=${1}
SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/functions.sh
source ${SCRIPT_DIR}/defaults.sh
dir=${dir_restart_ocean}
echo "DOWNLOADING MOM6 restarts to ${dir}"
mkdir -p ${dir} && cd ${dir}
aws_restart_path="${aws_path}/${dtg}/gdas.${dtg:0:8}/${dtg:8:10}/model/ocean/restart"
aws_inc_path="${aws_path}/${dtg}/gdas.${dtg:0:8}/${dtg:8:10}/analysis/ocean"

############
# MOM restarts
file_in=${aws_restart_path}/${DTG_TEXT_SRC}.MOM.res.nc  
WGET_AWS ${file_in} ${DTG_TEXT_DES}.MOM.res.nc
for i in $(seq 1 3); do
    file_in=${aws_restart_path}/${DTG_TEXT_SRC}.MOM.res_${i}.nc  
    WGET_AWS ${file_in} ${DTG_TEXT_DES}.MOM.res_${i}.nc
done
FIND_EMPTY_FILES ${PWD}

############
# MOM inc files
dir=${dir_inc_ocean}
mkdir -p ${dir} && cd ${dir}
echo "DOWNLOADING MOM6 increments to ${dir}"
file_in=${aws_inc_path}/gdas.t${dtg:8:10}z.ocn.incr.nc  
file_out=gdas.t${dtg:8:10}z.ocn.incr.nc
WGET_AWS ${file_in} ${file_out}
FIND_EMPTY_FILES ${PWD}

exit 0
