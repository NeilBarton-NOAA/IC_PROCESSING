#!/bin/bash
set -u
# Get ICs for GFS or CPC
dtg=${1:-2025070100}
export IC_SRC=REPLAY
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

####################################
files="ATM OCN"

####################################
BACKGROUND_JOB=T && export MV_DATA=T && export DOWNLOAD=T
for f in ${files}; do
    JOB_NAME=GETPERTS.${f}.${dtg} && echo ${JOB_NAME}
    source ${TOPDIR}/MACHINE/config.sh
    ${SUBMIT_HPSS} ${SCRIPT_DIR}/get_perturbations_${f}.sh ${dtg}  
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
    exit 1
done 

