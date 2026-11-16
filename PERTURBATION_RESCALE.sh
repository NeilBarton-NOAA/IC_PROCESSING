#!/bin/bash
set -u
export dtg=${1:-2024060100}
export IC_SRC=${2:-CPC_landice}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh
export OCEAN_PERT_ORIGIN=${IC_DIR}
JOB_NAME=PERT.RESCALE.${dtg}
WALLTIME="00:30:00"
BACKGROUND_JOB=T
source ${TOPDIR}/MACHINE/config.sh
echo "${JOB_NAME}"
if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
    ${SUBMIT} ${SCRIPT_DIR}/rescale_perts_OCN.sh ${dtg}
elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
    echo -e "${SUBMIT}\n${SCRIPT_DIR}/rescale_perts_OCN.sh ${dtg}" > submit_PERT_RESCALE.sh
    qsub submit_PERT_RESCALE.sh
    rm submit_PERT_RESCALE.sh
fi
[[ ${?} > 0 ]] && echo "FATAL with SUBMIT PERTURBATION RESCALE" && exit 1
