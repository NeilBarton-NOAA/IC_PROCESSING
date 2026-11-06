#!/bin/bash
set -u
export dtg=${1:-2024060100}
export IC_SRC=${2:-CPC_landice}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

members="000"
for mem in ${members}; do
    JOB_NAME=ICE.MEM${mem}.${dtg}
    WALLTIME="00:30:00"
    BACKGROUND_JOB=T
    source ${TOPDIR}/MACHINE/config.sh
    echo "${JOB_NAME}"
    if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
        ${SUBMIT} ${SCRIPT_DIR}/replace_ICE.sh ${dtg}
    elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
        echo -e "${SUBMIT}\n${SCRIPT_DIR}/replace_ICE.sh ${dtg}" > submit_ICE_REPLACE.sh
        qsub submit_ICE_REPLACE.sh
        rm submit_ICE_REPLACE.sh
    fi
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_ICE_REPLACE" && exit 1
done
