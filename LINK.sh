#!/bin/bash
set -ux
# Link ICs for GFS_DA ICs,CPC ICs (CPC_land, CPC_landice, CPC_landice) and Replay ICs
dtg=${1:-2022030100}
export IC_SRC=${2:-REPLAY}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh
BACKGROUND_JOB=F
JOB_NAME=LINK.${IC_SRC}.${dtg} && echo ${JOB_NAME}

if [[ ${IC_SRC} == "GFS_DA" ]]; then
# Remove unnecessary folders & files
    rm -rf ${IC_DIR}/${run}.${dtg:0:8}/${dtg:8:2}/mem000/model/atmos/restart
    rm -rf ${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/atmos
    rm -rf ${dir_inc_atmos} ${dir_inc_ocean}
fi

source ${TOPDIR}/MACHINE/config.sh
echo "${JOB_NAME}"
    if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
        ${SUBMIT} ${SCRIPT_DIR}/link_member_dirs.sh ${dtg} ${IC_SRC}
    elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
        echo -e "${SUBMIT}\n${SCRIPT_DIR}/link_member_dirs.sh ${dtg} ${IC_SRC}" > submit_LINK.sh
        qsub submit_LINK.sh
        rm submit_LINK.sh
    fi
[[ ${?} > 0 ]] && echo "FATAL with SUBMIT_LINK" && exit 1

echo "LINK ${IC_SRC} ICs is Complete"
