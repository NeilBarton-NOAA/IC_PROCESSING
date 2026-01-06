#!/bin/bash
set -u
dtg=${1}
SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/defaults.sh
inc_dir=${dir_inc_atmos}/perturbation
inc_dir=${inc_dir/mem000/mem001}
mkdir -p ${inc_dir} && cd ${inc_dir}

echo "DOWNLOADING ATM IC PERTURBATIONS TO ${inc_dir}"
if (( ${dtg:0:4} > 2024 )); then
    year=$(( ${dtg:0:4} - 10 ))
    fdtg=${year}${dtg:4:8}
    echo "Year after before 2025, grabbing 10 years before ${fdtg}"
else
    fdtg=${dtg}
fi

########################
# Atmosphere perturbation files on hpss
if [[ ${ATMRES} == "C96" ]]; then
    hpss_atm_increment_dir=/ESRL/BMC/gsienkf/Permanent/UFS_replay_input/era5/C96_perts
    if [[ ${dtg:3:1} == 0 ]]; then
        EY=1
    else
        EY=0
    fi
    file_name=${hpss_atm_increment_dir}/atm_perts_for_SFS_${ATMRES}_${EY}${dtg:3:3}01.tar
else # C192 or C384
    hpss_atm_increment_dir=/ESRL/BMC/gsienkf/Permanent/UFS_replay_input/era5/C384_perts/${dtg:0:4}
    if [[ "${ATMRES}" == "C192" ]]; then
        file_name=$( hsi -q ls -l ${hpss_atm_increment_dir}/C384_era5anl_${fdtg:0:6}*03_inc.tar 2>&1 | grep C384_era5anl | head -n 1 | awk '{print $9}' )
        file_name=${hpss_atm_increment_dir}/${file_name}
    else
        if [[ ${NENS} == 10 ]]; then
            file_name=${hpss_atm_increment_dir}/${ATMRES}_era5anl_${dtg:0:8}03_inc.tar
        else
            file_name=${hpss_atm_increment_dir}/${ATMRES}_era5anl_5mem_${dtg:0:8}03_inc.tar
        fi
    fi
fi

htar -xvf ${file_name}
if (( ${?} > 0 )); then
    echo 'ERROR in htar, file likely does not exist'
    echo '  file_name:', ${file_name}
    hpss_atm_increment_dir=/ESRL/BMC/gsienkf/2year/whitaker/era5/C384ensperts
    if [[ ${ATMRES} == "C192" ]]; then
        file_name=$( hsi -q ls -l ${hpss_atm_increment_dir}/C384_era5anl_${fdtg:0:6}*03_inc.tar 2>&1 | grep C384_era5anl | head -n 1 | awk '{print $9}' )
        NENS=10
        file_name=${hpss_atm_increment_dir}/${file_name}
    else
        file_name=${hpss_atm_increment_dir}/$(basename ${file_name})
    fi
    htar -xvf ${file_name}
    if (( ${?} > 0 )); then
        echo 'FATAL in htar, file also not at'
        echo '  file_name:', ${file_name}
        exit 1
    fi
fi

########################
# copy increment files to directories
orig_dir=${dir_inc_atmos/mem000/mem001}
for n in $( seq 1 ${NENS}); do
    # copy file to correct directory
    mem=$(printf "%03d" ${n})
    dir_mem=${orig_dir/mem001/mem${mem}}
    mkdir -p ${dir_mem}
    if [[ ${NENS} == 10 ]]; then
        i=$(( n - 1 ))
    else
        i=$(( n + 4 ))
    fi
    inc_file=${dir_mem}/${DTG_TEXT_DES}.fv3_perturbation.nc
    if [[ ${ATMRES} == "C96" ]]; then
        hpss_file=${inc_dir}/${EY}${dtg:3:3}01/${ATMRES}_era5anl_mem${mem}_${EY}${dtg:3:3}01.nc 
    else
        hpss_file=$( ls ${inc_dir}/C384_era5anl_???${i}_${fdtg:0:6}??03.nc | head -n 1 )
    fi
    mv ${hpss_file} ${inc_file}
    if (( ${?} > 0 )); then
        echo 'ERROR in copying perturbation'
        echo "  mv ${hpss_file} ${inc_file}"
        exit 1
    fi
done
rm -r ${inc_dir}
echo 'ATM IC PERTURBATION FILES DOWNLOADED AND PUT INTO MEM DIRECTORIES'
exit 0
