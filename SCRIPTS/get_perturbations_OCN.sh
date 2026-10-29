#!/bin/bash
#############
set -ux
dtg=${1}
IC_SRC=${2}
source "${SCRIPT_DIR}/defaults.sh"
source "${SCRIPT_DIR}/functions.sh"
inc_dir="${dir_inc_ocean}/perturbation"
mkdir -p "${inc_dir}" && cd "${inc_dir}"
echo "DOWNLOADING OCN IC PERTURBATION to ${inc_dir}"

##################################################################
## Download OCN perts on hpss

if (( ${dtg:0:4} > 2025 )); then
    year=$(( ${dtg:0:4} - 10 ))
    fdtg=${year}${dtg:4:8}
    echo "Year after before 2025, grabbing 10 years before ${fdtg}"
else
    fdtg=${dtg}
fi

# Ocean perturbation files for year ${fdtg:0:4} on hpss
hpss_ocn_increment_dir=/ESRL/BMC/gsienkf/Permanent/UFS_replay_input/oras5_ocn/ensemble_perts/${OCNRES}/${fdtg:0:4}
# Target date from your dtg variable (YYYYMMDDHH)
target_dtg="${fdtg}"
found=0
# Loop order: Exact date, then +1 day, then -1 day, +2, -2, etc.
# This ensures that if both +1 and -1 exist, the future (+1) is chosen first.
for offset in 0 1 -1 2 -2 3 -3 4 -4; do
    # Calculate the calendar date for the current offset
    current_date=$(date -d "${target_dtg:0:8} ${offset} days" +"%Y%m%d")
    potential_file="${hpss_ocn_increment_dir}/ocn_perts_10mem_${current_date}00.tar"

    # Check if file exists on HPSS
    if hsi -q "ls ${potential_file}"; then
        file_name="${potential_file}"
        found=1
        break
    fi
done

if [[ $found -eq 1 ]]; then
    echo "Selected file (Future Priority): ${file_name}"
else
    echo "No file found within +/- 4 days of ${target_dtg}"
    exit 1
fi

echo "DOWNLOADING: ${file_name}"
htar -xvf "${file_name}"
if (( ${?} > 0 )); then
    echo 'ERROR in htar, file likely does not exist'
    echo '  file_name:', ${file_name}
    exit 1
fi

#######################################################################
# Copy increment files to directories
orig_dir="${dir_inc_ocean/mem000/mem001}"
for n in $( seq 1 ${NENS} ); do
    # copy file to correct directory
    mem=$(printf "%03d" ${n})
    dir_mem="${orig_dir/mem001/mem${mem}}"
    mkdir -p "${dir_mem}"
    inc_file="${dir_mem}/${DTG_TEXT}.mom6_perturbation.nc"
    hpss_file="$( ls ${inc_dir}/??????????/mem${mem}_pert.nc )"
    echo "mv ${hpss_file} ${inc_file}"
    mv "${hpss_file}" "${inc_file}"
    if (( ${?} > 0 )); then
        echo 'ERROR in copying perturbation'
        echo "  mv ${hpss_file} ${inc_file}"
        exit 1
    fi
done

rm -r "${inc_dir}"
echo 'OCN IC PERTURBATION FILES DOWLOANDED AND PUT INTO MEM DIRS'

#########################################################################
## Download OCN perts from AWS

#LN=${NENS}
#aws_ocn_inc_dir="https://noaa-oar-sfsdev-pds.s3.amazonaws.com/index.html#input/ocn_ice/${OCNRES}/ens_perts"
#if (( ${dtg} > 2023110100 )); then
#    year=$(( ${dtg:0:4} - 10 ))
#    fdtg=${year}${dtg:4:8}
#    echo "Year is after 2024, grabbing 10 years before ${fdtg}"
#    else
#    fdtg=${dtg:0:8}
#fi

#for m in $(seq 1 10); do
#    mem=$(printf "%03d" ${m})
#    inc_dir=${dir_inc_ocean/mem000/mem${mem}}
#    mkdir -p ${inc_dir}
#    file_in=mem${mem}_pert.nc
#    file_out=${inc_dir}/${DTG_TEXT}.mom6_perturbation.nc
#    WGET_AWS ${aws_ocn_inc_dir}/${fdtg}/${file_in} ${file_out} 
#    [[ $? > 0 ]] && echo "FATAL in download" && exit 1
#done

#echo 'OCN IC PERTURBATION FILES DOWLOANDED AND PUT INTO MEM DIRS'
#################################################################################
