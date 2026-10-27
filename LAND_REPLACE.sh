#!/bin/bash
set -u
# Run REPLACE LAND STATES for ICs
export dtg=${1:-2022030100}
export IC_SRC=${2:-CPC_land}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh
file_land_states="${dir_land_states}/ufs_land_restart.${dtg:0:4}-${dtg:4:2}-${dtg:6:2}_00-00-00.nc"
PYTHON=$(which python)

members="000"
for mem in ${members}; do
    JOB_NAME=LAND.MEM${mem}.${dtg}
    WALLTIME="00:30:00"
    BACKGROUND_JOB=F
    source ${TOPDIR}/MACHINE/config.sh
    echo "${JOB_NAME}"
    if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
        ${SUBMIT} --wrap="${PYTHON} ${SCRIPT_DIR}/LAND_UPDATE/replace_land_states.py \"${dir_input_atmos}/\" ${file_land_states}" 
    elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
        echo -e "${SUBMIT}\n${PYTHON} ${SCRIPT_DIR}/LAND_UPDATE/replace_land_states.py \"${dir_input_atmos}/\" ${file_land_states}" > submit_LAND_REPLACE.sh
        qsub submit_LAND_REPLACE.sh
        rm submit_LAND_REPLACE.sh
    fi
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_LAND_REPLACE" && exit 1
done
