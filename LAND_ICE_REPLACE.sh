#!/bin/bash
set -u
# Run REPLACE LAND AND ICE STATES for ICs
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
    BACKGROUND_JOB=F
    source ${TOPDIR}/MACHINE/config.sh
    ## Conduct both land and ice state replacement for CPC_landice
    if [[ ${IC_SRC} == "CPC_landice" ]]; then
        JOB_NAME=LAND.MEM${mem}.${dtg}
        WALLTIME="00:30:00"
        echo "${JOB_NAME}"
        if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
            ${SUBMIT} --wrap="${PYTHON} ${SCRIPT_DIR}/LAND_UPDATE/replace_land_states.py \"${dir_input_atmos}/\" ${file_land_states}" 
        elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
            echo -e "${SUBMIT}\n${PYTHON} ${SCRIPT_DIR}/LAND_UPDATE/replace_land_states.py \"${dir_input_atmos}/\" ${file_land_states}" > submit_LAND_REPLACE.sh
            qsub submit_LAND_REPLACE.sh
            rm submit_LAND_REPLACE.sh
        fi
        [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_LAND_REPLACE" && exit 1

        JOB_NAME=ICE.MEM${mem}.${dtg}
        WALLTIME="00:30:00"
        echo "${JOB_NAME}"
        if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
            ${SUBMIT} ${SCRIPT_DIR}/ICE_UPDATE/run_cice6_insertion.sh --rdates ${dtg:0:8}
        elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
            echo -e "${SUBMIT}\n${SCRIPT_DIR}/ICE_UPDATE/run_cice6_insertion.sh --rdates ${dtg:0:8}" > submit_ICE_REPLACE.sh
            qsub submit_ICE_REPLACE.sh
            rm submit_ICE_REPLACE.sh
        fi
        [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_ICE_REPLACE" && exit 1
    else

    ## Only conduct land state replacement
        JOB_NAME=LAND.MEM${mem}.${dtg}
        WALLTIME="00:30:00"
        echo "${JOB_NAME}"
        if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
            ${SUBMIT} --wrap="${PYTHON} ${SCRIPT_DIR}/LAND_UPDATE/replace_land_states.py \"${dir_input_atmos}/\" ${file_land_states}"
        elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
            echo -e "${SUBMIT}\n${PYTHON} ${SCRIPT_DIR}/LAND_UPDATE/replace_land_states.py \"${dir_input_atmos}/\" ${file_land_states}" > submit_LAND_REPLACE.sh
            qsub submit_LAND_REPLACE.sh
            rm submit_LAND_REPLACE.sh
        fi
        [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_LAND_REPLACE" && exit 
    fi
done
