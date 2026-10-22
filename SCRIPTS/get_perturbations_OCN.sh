#!/bin/bash
#############
# HPSS location
#   hpss_ocn_increment_dir=/ESRL/BMC/gsienkf/Permanent/UFS_replay_input/oras5_ocn/ensemble_perts/${OCNRES}
#   file_name=${hpss_ocn_increment_dir}/ocn_perts_for_SFS_${OCNRES}_${dtg:0:6}0100.tar
set -ux
dtg=${1}
IC_SRC=${2}
source ${SCRIPT_DIR}/defaults.sh
source ${SCRIPT_DIR}/functions.sh
echo "DOWNLOADING OCN IC PERTURBATION"

############
# Ocean perturbation files on hpss
LN=${NENS}
aws_ocn_inc_dir="https://noaa-oar-sfsdev-pds.s3.amazonaws.com/index.html#input/ocn_ice/${OCNRES}/ens_perts"
if (( ${dtg} > 2023110100 )); then
    year=$(( ${dtg:0:4} - 10 ))
    fdtg=${year}${dtg:4:8}
    echo "Year is after 2024, grabbing 10 years before ${fdtg}"
    else
    fdtg=${dtg:0:8}
fi

for m in $(seq 1 10); do
    mem=$(printf "%03d" ${m})
    inc_dir=${dir_inc_ocean/mem000/mem${mem}}
    mkdir -p ${inc_dir}
    file_in=mem${mem}_pert.nc
    file_out=${inc_dir}/${DTG_TEXT}.mom6_perturbation.nc
    WGET_AWS ${aws_ocn_inc_dir}/${fdtg}/${file_in} ${file_out} 
    [[ $? > 0 ]] && echo "FATAL in download" && exit 1
done

echo 'OCN IC PERTURBATION FILES DOWLOANDED AND PUT INTO MEM DIRS'
