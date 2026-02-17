#!/bin/bash
set -u
# Get ICs for GFS or CPC
dtg=${1:-2025070100}
export IC_SRC=GFS 
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh
if [[ ${IC_SRC} == "GFS" ]]; then
    files="gdasocean_restart gdasocean_analysis gdas_restarta gdas_restartb enkfgdas_restarta_grp1 enkfgdas_restartb_grp1"
    files="gdasocean_restart gdasocean_analysis gdas_restarta gdas_restartb"
    files="enkfgdas_restartb_grp1"
    BACKGROUND_JOB=T && export MV_DATA=F && export DOWNLOAD=F
    for f in ${files}; do
       JOB_NAME=GET.${f}.${dtg} && echo ${JOB_NAME}
       source ${TOPDIR}/MACHINE/config.sh
       ${SUBMIT_HPSS} ${SCRIPT_DIR}/get_GFS.sh ${dtg} ${f} 
       [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
    done 
fi

#if [[ ${IC_SRC} == "GFS_EXACT_DATE" ]]
#
#fi
