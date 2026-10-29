#!/bin/bash
set -u
# Get ICs for GFS or CPC
dtg=${1:-2026030100}
export IC_SRC=${2:-GFS}
export ATMRES=${3:-C192}
export TOPDIR=${PWD}
export SCRIPT_DIR=${TOPDIR}/SCRIPTS
source ${TOPDIR}/MACHINE/config.sh
source ${TOPDIR}/SCRIPTS/defaults.sh
BACKGROUND_JOB=F && export MV_DATA=T && export DOWNLOAD=T

####################################
if [[ ${IC_SRC} == "GFS" ]]; then
    files="gdasocean_restart gdasocean_analysis gdas_restarta gdas_restartb enkfgdas \
           enkfgdas_restarta_grp1 enkfgdas_restartb_grp1"
    # comment out below if only need 10 members
    files="${files} enkfgdas_restarta_grp2 enkfgdas_restartb_grp2 enkfgdas_restarta_grp3 enkfgdas_restartb_grp3"
elif [[ ${IC_SRC} == "GFS_3DVAR" ]]; then
    files="gdasocean_restart gdasocean_analysis gdas_restarta gdas_restartb"
elif [[ ${IC_SRC} == "REPLAY" ]] || [[ ${IC_SRC} == *"CPC"* ]]; then
    files="ATM OCN ICE"
    [[ ${ATMRES} == "C384" ]] && files="${files} MED"
fi
[[ ${IC_SRC} == *"CPC"* ]] && SCRIPT_TAG=CPC
[[ ${IC_SRC} == "GFS_3DVAR"* ]] && SCRIPT_TAG=GFS
####################################
for f in ${files}; do
    JOB_NAME=GET.${IC_SRC}.${f}.${dtg} && echo ${JOB_NAME}
    source ${TOPDIR}/MACHINE/config.sh
    if [[ ${BATCH_SYSTEM} == "sbatch" ]] || [[ ${BACKGROUND_JOB} == T ]]; then
        ${SUBMIT_HPSS} ${SCRIPT_DIR}/get_${SCRIPT_TAG:-$IC_SRC}.sh ${dtg} ${f} ${IC_SRC} 
    elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
        echo -e "${SUBMIT_HPSS}\n${SCRIPT_DIR}/get_${SCRIPT_TAG:-$IC_SRC}.sh ${dtg} ${f} ${IC_SRC}" > submit_GETICS.sh
        qsub submit_GETICS.sh
        rm submit_GETICS.sh
    fi
    [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS_GETICS" && exit 1
done 

# GET REQUIRED PERTURBATIONS
####################################
if [[ ${IC_SRC} == "CPC_land" || ${IC_SRC} == "CPC_landice" || ${IC_SRC} == "GFS_3DVAR" || ${IC_SRC} == "REPLAY" ]]; then
    files="ATM OCN"
    for f in ${files}; do
        JOB_NAME=GETPERTS.${IC_SRC}.${f}.${dtg} && echo ${JOB_NAME}
        source ${TOPDIR}/MACHINE/config.sh
        ${SUBMIT_HPSS} ${SCRIPT_DIR}/get_perturbations_${f}.sh ${dtg} ${IC_SRC}
        [[ ${?} > 0 ]] && echo "FATAL with SUBMIT_HPSS_GETPERTS" && exit 1
    done
fi
