#!/bin/bash
set -u
dtg=${1}
SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/defaults.sh
inc_dir=${dir_inc_atmos}/perturbation
inc_dir=${inc_dir/mem000/mem001}
mkdir -p ${inc_dir} && cd ${inc_dir}

echo "DOWNLOADING ATM IC PERTURBATIONS TO ${inc_dir}"
########################
# Atmosphere perturbation files on hpss
if [[ ${ATMRES} == "C96" ]]; then
    hpss_atm_increment_dir=/ESRL/BMC/gsienkf/Permanent/UFS_replay_input/era5/C96_perts
    [[ ${dtg:3:1} == 0 ]] && EY=1
    file_name=${hpss_atm_increment_dir}/atm_perts_for_SFS_${ATMRES}_${EY:-0}${dtg:3:3}01.tar
else # C192 or C384
    hpss_path1=/ESRL/BMC/gsienkf/Permanent/UFS_replay_input/era5/C384_perts
    hpss_path2=/ESRL/BMC/gsienkf/2year/whitaker/era5/C384ensperts
    files=$( hsi find ${hpss_path1} ${hpss_path2} -name "C384_era5anl_[12]???*${dtg:4:4}??_inc.tar" 2>&1 | grep tar )
    years=$( grep -oP '/\d{4}/' <<< "${files}" | tr -d '/' | sort -u ) 
    if echo ${years} | grep -q ${dtg:0:4}; then
        file_name=$( grep ${dtg:0:8} <<< ${files} )
    else
        echo "WARNING year not found"
        file_name=$( head -n 1 <<< ${files} )
    fi
fi

echo "htar -xvf ${file_name}"
#htar -xvf ${file_name}
if (( ${?} > 0 )); then
    echo 'ERROR in htar, file likely does not exist'
    echo '  file_name:', ${file_name}
fi
########################
# copy increment files to directories
orig_dir=${dir_inc_atmos/mem000/mem001}
for n in $( seq 1 ${NENS} ); do
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
        hpss_file=$( ls ${inc_dir}/C384_era5anl_inc${i}_*.nc )
    fi
    echo "mv ${hpss_file} ${inc_file}"
    mv ${hpss_file} ${inc_file}
    if (( ${?} > 0 )); then
        echo 'ERROR in copying perturbation'
        echo "  mv ${hpss_file} ${inc_file}"
        exit 1
    fi
done
rm -r ${inc_dir}
echo 'SUCCESFUL: ATM IC PERTURBATION FILES DOWNLOADED AND PUT INTO MEM DIRECTORIES'
exit 0
