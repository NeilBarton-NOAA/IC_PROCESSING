#!/bin/bash
set -u
# Get ICs for GFS or CPC
dtg=${1:-2022030100}
export IC_SRC=${2:-CPC_land} #CPC_land, CPC_landice,GFS_3DVAR,REPLAY
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

####################################
files="ATM OCN"

####################################
BACKGROUND_JOB=F
for f in ${files}; do
    JOB_NAME=GETPERTS.${IC_SRC}.${f}.${dtg} && echo ${JOB_NAME}
    source ${TOPDIR}/MACHINE/config.sh
    ${SUBMIT_HPSS} ${SCRIPT_DIR}/get_perturbations_${f}.sh ${dtg} ${IC_SRC}
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
done 

