#!/bin/bash
set -u
# Link ICs for GFS_DA ICs,CPC ICs (CPC_land, CPC_landice, CPC_landice) and Replay ICs
dtg=${1:-2022030100}
export IC_SRC=${2:-REPLAY}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

${SCRIPT_DIR}/link_member_dirs.sh ${dtg} ${IC_SRC} 

if [[ ${IC_SRC} == "GFS_DA" ]]; then
# Remove unnecessary folders & files
rm -rf ${dir_restart_atmos} ${dir_inc_atmos} ${dir_inc_ocean}
echo "Finalizing GFS_DA ICs Dir is Complete"
fi

