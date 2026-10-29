#!/bin/bash
set -u
export dtg=${1:-2022030100}
export IC_SRC=${2:-CPC_land}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh
export rdate_ice=${dtg_minus6:0:8}
export hr_ice=${dtg_minus6:8:2}

members="000"
for mem in ${members}; do
    JOB_NAME=ICE.MEM${mem}.${dtg}
    WALLTIME="00:30:00"
    BACKGROUND_JOB=F
    source ${TOPDIR}/MACHINE/config.sh
    echo "${JOB_NAME}"
    if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
        ${SUBMIT} ${SCRIPT_DIR}/ICE_UPDATE/run_cice6_insertion.sh --rdates ${dtg:0:8}
    elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
        echo -e "${SUBMIT}\n${SCRIPT_DIR}/ICE_UPDATE/run_cice6_insertion.sh --rdates ${dtg:0:8}" > submit_ICE_REPLACE.sh
        qsub submit_ICE_REPLACE.sh
        rm submit_ICE_REPLACE.sh
    fi
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_ICE_REPLACE" && exit 1
done
