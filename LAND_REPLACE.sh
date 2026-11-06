#!/bin/bash
set -u
# Run REPLACE LAND STATES for ICs
export dtg=${1:-2022030100}
export IC_SRC=${2:-CPC_land}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

members="000"
for mem in ${members}; do
    JOB_NAME=LAND.MEM${mem}.${dtg}
    WALLTIME="00:30:00"
    BACKGROUND_JOB=F
    source ${TOPDIR}/MACHINE/config.sh
    echo "${JOB_NAME}"
    if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
        ${SUBMIT} ${SCRIPT_DIR}/replace_LND.sh ${dtg}
    elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
        echo -e "${SUBMIT}\n${SCRIPT_DIR}/replace_LND.sh ${dtg}" > submit_LAND_REPLACE.sh
        qsub submit_LAND_REPLACE.sh
        rm submit_LAND_REPLACE.sh
    fi
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_LAND_REPLACE" && exit 1
done
