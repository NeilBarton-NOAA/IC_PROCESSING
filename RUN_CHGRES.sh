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

set -x
if [[ ${IC_SRC} == "GFS" ]]; then
    models="ATM"
    members=$( ls -d ${dir_restart_atmos%/mem000*}/mem*/ | grep -oP '(?<=mem)\d{3}')
fi
exit 1
for model in ${models}; do
    for mem in ${members}; do
        echo ${mem}
        JOB_NAME=CHGRES.${model}.MEM${mem}.${dtg}
        NTASKS=12
        WALLTIME="00:30:00"
        source ${TOPDIR}/MACHINE/config.sh
        if [[ ${mem} == "000" ]]; then
            ATMRES="C1152"
        else
            ATMRES="C384"
        fi
        ${SUBMIT} ${SCRIPT_DIR}/chgres_${model}.sh ${dtg} ${ATMRES} mx025 ${mem}
        [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
    done
done
