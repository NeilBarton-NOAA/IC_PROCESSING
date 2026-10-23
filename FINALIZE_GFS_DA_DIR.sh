#!/bin/bash
set -ux
dtg=${1:-2026030100}
export IC_SRC=${2:-GFS_DA}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

####################################
# Remove unnecessary folders & files
rm -rf ${dir_restart_atmos} ${dir_inc_atmos} ${dir_inc_ocean}
echo "Finalizing GFS_DA ICs Dir is Complete"
