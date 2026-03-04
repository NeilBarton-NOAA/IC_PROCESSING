#!/bin/bash
set -u
# Get ICs for GFS or CPC
dtg=${1:-2026030100}
export IC_SRC=GFS 
#export IC_SRC=REPLAY
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh

####################################
if [[ ${IC_SRC} == "GFS" ]]; then
    files="gdasocean_restart gdasocean_analysis gdas_restarta gdas_restartb enkfgdas \
           enkfgdas_restarta_grp1 enkfgdas_restartb_grp1 \
           enkfgdas_restarta_grp2 enkfgdas_restartb_grp2 \
           enkfgdas_restarta_grp3 enkfgdas_restartb_grp3"
elif [[ ${IC_SRC} == "REPLAY" ]]; then
    files="ATM OCN ICE"
    files="ICE"
fi

####################################
BACKGROUND_JOB=F && export MV_DATA=T && export DOWNLOAD=T
for f in ${files}; do
    JOB_NAME=GET.${f}.${dtg} && echo ${JOB_NAME}
    source ${TOPDIR}/MACHINE/config.sh
    ${SUBMIT_HPSS} ${SCRIPT_DIR}/get_${IC_SRC}.sh ${dtg} ${f} 
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS" && exit 1
done 

