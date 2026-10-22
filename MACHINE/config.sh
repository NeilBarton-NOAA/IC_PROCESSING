#!/bin/sh
####################################
# Local Directoriesi and 'global' variables
export HPC_ACCOUNT=${HPC_ACCOUNT:-"marine-cpu"}
export TOPDIR=${SCRIPT_DIR}/../
machine=$(uname -n)
JOB_NAME=${JOB_NAME:-hpss}
WALLTIME=${WALLTIME:-01:00:00}
NTASKS=${NTASKS:-1}

BATCH_SYSTEM="sbatch"
SUBMIT_SUFFIX=""
SUBMIT_HPSS_SUFFIX=""
export APRUN="srun"
############
# scrath dir based in machine
machine=$(uname -n)
export chgres_compiler=intel
if [[ ${machine} == hercules* ]]; then
    machine=hercules
    export WORK_DIR=/work/noaa/marine/${USER}
    CODE_DIR=${WORK_DIR}/CODE/IC_PROCESSING
    export CICE6_ICinsertion_DATA=${CODE_DIR}/CICE6_ICinsertion_DATA
elif [[ ${machine} == gaea* || ${machine} == dtn* || ${machine} == c6* ]]; then
    machine=gaea
    export WORK_DIR=/gpfs/f6/sfs-emc/scratch/${USER}
    CODE_DIR=${WORK_DIR}/CODE/IC_PROCESSING
    export CICE6_ICinsertion_DATA=${CODE_DIR}/CICE6_ICinsertion_DATA
    SUBMIT_SUFFIX="--qos=normal --clusters=c6 --partition=batch"
    SUBMIT_HPSS_SUFFIX="--mem=100G --qos=hpss --clusters=es --partition=dtn_f5_f6 --constraint=f6"
elif [[ ${machine} == u* ]]; then
    machine=ursa
    export chgres_compiler=intelllvm
    export WORK_DIR=/scratch4/NCEPDEV/stmp/${USER}
    CODE_DIR=${WORK_DIR}/CODE/IC_PROCESSING
    export CICE6_ICinsertion_DATA=${CODE_DIR}/CICE6_ICinsertion_DATA
    SUBMIT_SUFFIX="--mem=0 --qos=batch"
    SUBMIT_HPSS_SUFFIX="--mem=100G --partition=u1-service"
elif [[ ${machine} == *[cd]login* ]] || [[ ${machine} == nid* ]] || [[ ${machine} == *dx* ]]; then
    BATCH_SYSTEM="qsub"
    machine=wcoss2
    CODE_DIR=/lfs/h2/emc/couple/noscrub/${USER}/CODE/IC_PROCESSING
    export CICE6_ICinsertion_DATA=${CODE_DIR}/CICE6_ICinsertion_DATA
    export WORK_DIR=/lfs/h2/emc/stmp/${USER}
    export APRUN="mpiexec"
else
    echo 'FATAL: MACHINE UNKNOWN', ${machine}
    exit 1
fi
    
if [[ ${BATCH_SYSTEM} == "sbatch" ]]; then
SUBMIT_BASE="${BATCH_SYSTEM} 
    --job-name=${JOB_NAME} 
    --output=${TOPDIR}/logs/${JOB_NAME}.out
    --error=${TOPDIR}/logs/${JOB_NAME}.out
    --time=${WALLTIME} 
    --account=${HPC_ACCOUNT} 
    --ntasks=${NTASKS}"

SUBMIT="${SUBMIT_BASE} ${SUBMIT_SUFFIX}"
SUBMIT_HPSS="${SUBMIT_BASE} ${SUBMIT_HPSS_SUFFIX}"

elif [[ ${BATCH_SYSTEM} == "qsub" ]]; then
SUBMIT="#!/bin/bash
#PBS -N ${JOB_NAME} 
#PBS -o $(readlink -m ${TOPDIR}/logs/${JOB_NAME}.out) 
#PBS -e $(readlink -m ${TOPDIR}/logs/${JOB_NAME}.out) 
#PBS -l walltime=0${WALLTIME}
#PBS -l select=1:ncpus=${NTASKS} 
#PBS -A ${HPC_ACCOUNT}
#PBS -V"
SUBMIT_HPSS="#!/bin/bash
#PBS -N ${JOB_NAME} 
#PBS -o $(readlink -m ${TOPDIR}/logs/${JOB_NAME}.out) 
#PBS -e $(readlink -m ${TOPDIR}/logs/${JOB_NAME}.out) 
#PBS -l walltime=0${WALLTIME}
#PBS -l select=1:ncpus=1 
#PBS -A ${HPC_ACCOUNT}
#PBS -q dev_transfer
#PBS -V"
fi

if [[ ${BACKGROUND_JOB:-F} == T ]]; then
    SUBMIT=""
    SUBMIT_HPSS=""
fi
mkdir -p ${TOPDIR}/logs
source ${SCRIPT_DIR}/../MACHINE/modules.sh
