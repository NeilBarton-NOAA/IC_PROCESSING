#!/bin/bash
set -u
# Run CHGRES for ICs
dtg=${1:-2025070100}
export DEBUG=F
export IC_SRC=GFS 
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

if [[ ${IC_SRC} == "GFS" ]]; then
    members=$( ls -d ${dir_inc_atmos%/mem000*}/mem*/ | grep -oP '(?<=mem)\d{3}')
fi
for mem in ${members}; do
    ATMRES="C384"
    JOB_NAME=REGRID.SFC.MEM${mem}.${dtg}
    NTASKS=12
    WALLTIME="00:30:00"
    BACKGROUND_JOB=F
    source ${TOPDIR}/MACHINE/config.sh
    echo "${JOB_NAME}"
    ${SUBMIT} ${SCRIPT_DIR}/regrid_SFC.sh ${dtg} ${ATMRES} ${mem}
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
done
