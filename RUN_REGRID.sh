#!/bin/bash
set -u
# Run REGRID for ICs
dtg=${1:-2026030100}
export IC_SRC=${2:-GFS} 
export DEBUG=F
export TOPDIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

if [[ ${IC_SRC} == "GFS" ]]; then
    members=$( ls -d ${dir_inc_atmos%/mem000*}/mem*/ | grep -oP '(?<=mem)\d{3}')
fi
for mem in ${members}; do
    ATMRES="C384"
    JOB_NAME=REGRID.${IC_SRC}.SFC.MEM${mem}.${dtg}
    NTASKS=12
    WALLTIME="02:00:00"
    BACKGROUND_JOB=F
    source ${TOPDIR}/MACHINE/config.sh
    echo "${JOB_NAME}"
    if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
        ${SUBMIT} ${SCRIPT_DIR}/regrid_SFC.sh ${dtg} ${ATMRES} ${mem}
    elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
        echo -e "${SUBMIT}\n${SCRIPT_DIR}/regrid_SFC.sh ${dtg} ${ATMRES} ${mem}" > submit_REGRID.sh
        qsub submit_REGRID.sh
        rm submit_REGRID.sh
    fi
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
done
