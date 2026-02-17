#!/bin/sh
####################################
# Local Directoriesi and 'global' variables
export HPC_ACCOUNT=${COMPUTE_ACCOUNT}
machine=$(uname -n)
JOB_NAME=${JOB_NAME:-hpss}
WALLTIME=${WALLTIME:-01:00:00}
NTASKS=${NTASKS:-1}

BATCH_SYSTEM="sbatch"
SUBMIT_SUFFIX=""
SUBMIT_HPSS_SUFFIX=""
if [[ ${machine:0:3} == hfe || ${machine} == h*[cm]* ]]; then
    machine=hera
    export WORK_DIR=/scratch2/NCEPDEV/stmp3
elif [[ ${machine} == hercules* ]]; then
    machine=hercules
    export WORK_DIR=/work/noaa/marine
elif [[ ${machine} == gaea* || ${machine} == dtn* || ${machine} == c6* ]]; then
    machine=gaea
    export WORK_DIR=/gpfs/f6/sfs-emc/scratch
    SUBMIT_SUFFIX="--qos=normal --clusters=c6 --partition=batch"
    SUBMIT_HPSS_SUFFIX="--mem=100G --qos=hpss --clusters=es --partition=dtn_f5_f6 --constraint=f6"
elif [[ ${machine} == u* ]]; then
    machine=ursa
    export WORK_DIR=/scratch4/NCEPDEV/stmp
    SUBMIT_SUFFIX="--mem=0 --qos=debug"
    SUBMIT_HPSS_SUFFIX="--mem=100G --partition=u1-service"
else
    echo 'FATAL: MACHINE UNKNOWN'
    exit 1
fi
SUBMIT_BASE="${BATCH_SYSTEM} 
    --job-name=${JOB_NAME} 
    --output=${TOPDIR}/logs/${JOB_NAME}.out
    --error=${TOPDIR}/logs/${JOB_NAME}.out
    --time=${WALLTIME} 
    --account=${HPC_ACCOUNT} 
    --ntasks=${NTASKS}"
SUBMIT="${SUBMIT_BASE} ${SUBMIT_SUFFIX}"
SUBMIT_HPSS="${SUBMIT_BASE} ${SUBMIT_HPSS_SUFFIX}"

if [[ ${BACKGROUND_JOB:-F} == T ]]; then
    SUBMIT=""
    SUBMIT_HPSS=""
fi

