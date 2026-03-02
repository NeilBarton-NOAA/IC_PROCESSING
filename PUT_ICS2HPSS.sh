#!/bin/bash
set -u
# Run CHGRES for ICs
dtg=${1:-2025070100}
export DEBUG=F
export IC_SRC=GFS 
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
BACKGROUND_JOB=F
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh
########################
JOB_NAME=PUT.HPSS.${dtg}
source ${TOPDIR}/MACHINE/config.sh
echo "${JOB_NAME}"
${SUBMIT_HPSS} ${SCRIPT_DIR}/hpss_put_ICs.sh ${dtg} 
[[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
