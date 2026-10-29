#!/bin/bash
set -ux
# Run CHGRES for ICs
dtg=${1:-2026030100}
export IC_SRC=${2:-GFS}
export DEBUG=F
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
BACKGROUND_JOB=F
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

if [[ ${IC_SRC} == "GFS" ]]; then
    models="ATM"
    members=$( ls ${dir_restart_atmos%/mem000*}/mem*/model/atmos/restart/*fv_core.res.nc | grep -oP '(?<=mem)\d{3}')
elif [[ ${IC_SRC} == "REPLAY" || ${IC_SRC} == "GFS_3DVAR" ]]; then
    models="ATM"
    members="000"
elif [[ ${IC_SRC} == *"CPC"* ]]; then
    models="ATM OCN ICE"
    members="000"
fi

for model in ${models}; do
    for mem in ${members}; do
        JOB_NAME=CHGRES.${IC_SRC}.${model}.MEM${mem}.${dtg}
        mem=${mem}
        if [[ ${mem} == "000" && ${IC_SRC} == GFS ]]; then
            ATMRES="C1152"
            WALLTIME="01:00:00"
            NTASKS=12 # WCOSS2
            #NTASKS=24
        else
            ATMRES="C384"
            WALLTIME="00:30:00"
            NTASKS=12
        fi
        source ${TOPDIR}/MACHINE/config.sh
        echo "${JOB_NAME}"
        if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
            ${SUBMIT} ${SCRIPT_DIR}/chgres_${model}.sh ${dtg} ${ATMRES} mx025 ${mem} ${IC_SRC} ${NTASKS}
        elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
            echo -e "${SUBMIT}\n${SCRIPT_DIR}/chgres_${model}.sh ${dtg} ${ATMRES} mx025 ${mem} ${IC_SRC} ${NTASKS}" > submit_CHGRES.sh
            qsub submit_CHGRES.sh
            rm submit_CHGRES.sh
        fi
        [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
    done
done
