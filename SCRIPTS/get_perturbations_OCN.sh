#!/bin/bash
set -u
dtg=${1}
SCRIPT_DIR=$(dirname "$0")
source ${SCRIPT_DIR}/defaults.sh
inc_dir=${dir_inc_ocean}/perturbation
inc_dir=${inc_dir/mem000/mem001}
mkdir -p ${inc_dir} && cd ${inc_dir}
echo "DOWNLOADING OCN IC PERTURBATION TO ${inc_dir}"

############
# Ocean perturbation files on hpss
LN=${NENS}
#aws_ocn_increment_dir="https://noaa-oar-sfsdev-pds.s3.amazonaws.com/input/ocn_ice/mx100/ens_perts"
hpss_ocn_increment_dir=/ESRL/BMC/gsienkf/Permanent/UFS_replay_input/oras5_ocn/ensemble_perts/${OCNRES}
if (( ${dtg:0:4} > 2023 )); then
    year=$(( ${dtg:0:4} - 10 ))
    dtg=${year}${dtg:4:8}
    echo "Year is after 2024, grabbing 10 years before ${dtg}"
fi
file_name=${hpss_ocn_increment_dir}/ocn_perts_for_SFS_${OCNRES}_${dtg:0:6}0100.tar
echo "DOWNLOADING: ${file_name}"
htar -xvf ${file_name}
if (( ${?} > 0 )); then
    echo 'ERROR in htar, file likely does not exist'
    echo '  file_name:', ${file_name}
    exit 1
fi

########################
# copy increment files to directories
orig_dir=${dir_inc_ocean/mem000/mem001}
for n in $( seq 1 ${LN} ); do
    # copy file to correct directory
    mem=$(printf "%03d" ${n})
    dir_mem=${orig_dir/mem001/mem${mem}}
    mkdir -p ${dir_mem}
    pert_file=${inc_dir}/??????????/mem${mem}_pert.nc
    inc_file=${dir_mem}/${DTG_TEXT_DES}.mom6_perturbation.nc
    mv ${pert_file} ${inc_file}
    if (( ${?} > 0 )); then
        echo 'ERROR in copying perturbation'
        echo "  mv ${pert_file} ${inc_file}"
        exit 1
    fi
done
exit 1
if [[ ${LN} == 4 ]]; then
    dir_mem=${dir_inc_ocean/mem001/mem005}
    dir_mem001=${dir_inc_ocean}
    mkdir -p ${dir_mem} && cd ${dir_mem}
    cp ${dir_mem001}/${DTG_TEXT_DES}.mom6_perturbation.nc .
fi
rm -r ${inc_dir}
echo 'OCN IC PERTURBATION FILES DOWLOANDED AND PUT INTO MEM DIRS'
