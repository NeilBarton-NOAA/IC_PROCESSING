#!/bin/bash
set -u
dtg=${1}
SCRIPT_DIR=$(dirname "$0")
SCOUT_RESTARTS=${SCOUT_RESTARTS:-T}
source ${SCRIPT_DIR}/functions.sh
source ${SCRIPT_DIR}/defaults.sh
restart_tile_files='ca_data fv_core.res fv_srf_wnd.res fv_tracer.res phy_data sfc_data'
restart_nontile_files='ca_data fv_core.res'
analysis_files="gdas.t${dtg:8:10}z.atmanl.nc gdas.t${dtg:8:10}z.sfcanl.nc"
dir=${dir_restart_atmos}
if [[ ${IC_SRC} == "SCOUT" ]]; then
    if [[ ${SCOUT_RESTARTS} == "T" ]]; then
        aws_path="${aws_path}/${dtg}/gdas.${dtg:0:8}/${dtg:8:10}/model/atmos/restart"
    else
        dir=${dir_input_atmos}
        aws_path="${aws_path}/${dtg}/gdas.${dtg:0:8}/${dtg:8:10}/analysis/atmos"
    fi
fi
mkdir -p ${dir} && cd ${dir}
echo "DOWNLOADING FV3 data to ${dir}"

if [[ ${IC_SRC} == "SCOUT" ]] && [[ ${SCOUT_RESTARTS} != "T" ]]; then
    for f in ${analysis_files}; do
        file_in=${aws_path}/${f}
        file_out=${f}
        WGET_AWS ${file_in} ${file_out} 
    done
else
    for f in ${restart_tile_files}; do
    for tile in $(seq 1 6); do
        file_in=${aws_path}/${DTG_TEXT_SRC}.${f}.tile${tile}.nc 
        file_out=${DTG_TEXT_DES}.${f}.tile${tile}.nc
        WGET_AWS ${file_in} ${file_out} 
    done
    done
    for f in ${restart_nontile_files}; do
        file_in=${aws_path}/${DTG_TEXT_SRC}.${f}.nc 
        file_out=${DTG_TEXT_DES}.${f}.nc
        WGET_AWS ${file_in} ${file_out} 
    done
    for tile in $(seq 1 6); do
        file_out=${DTG_TEXT_DES}.sfc_data.tile${tile}.nc
        ncatted -a checksum,,d,, ${file_out}
    done
    #if [[ ${ATMRES} == "C384" ]]; then
    #    touch ${DTG_TEXT_DES}.coupler.res
    #fi
fi

FIND_EMPTY_FILES ${PWD}

