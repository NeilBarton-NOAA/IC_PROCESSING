#!/bin/bash
set -u
export dtg=${1}
model=${2} # ATM, OCN, ICE
DEBUG=${DEBUG:-F}
export IC_SRC=REPLAY
source ${SCRIPT_DIR}/functions.sh
source ${SCRIPT_DIR}/defaults.sh

[[ ${model} == "ATM" ]] && dir=${dir_restart_atmos} && files=${restart_nontile_files_atmos}
[[ ${model} == "OCN" ]] && dir=${dir_restart_ocean} && files=${restart_files_ocean}
[[ ${model} == "ICE" ]] && dir=${dir_restart_ice} && files="iced.${dtg:0:4}-${dtg:4:2}-${dtg:6:2}-10800"

############
mkdir -p ${dir} && cd ${dir}
echo "DOWNLOADING ${model} data to ${dir}"

if [[ ${model} == "ATM" ]]; then
    for f in ${restart_tile_files_atmos}; do
    for tile in $(seq 1 6); do
        file_in=${f}.tile${tile}.nc 
        file_out=${DTG_TEXT}.${f}.tile${tile}.nc
        if [[ ${f} == "sfc_data" ]] && [[ "${dtg}" -ge "1994050100" && "${DTG}" -le "2023110106" ]]; then
            WGET_AWS ${aws_path_sfc}/${file_in} ${file_out} 
            [[ $? > 0 ]] && echo "FATAL in download" && exit 1
        else
            WGET_AWS ${aws_path_replay}/${file_in} ${file_out} 
            [[ $? > 0 ]] && echo "FATAL in download" && exit 1
        fi
    done
    done
fi
for f in ${files}; do
    file_in=${f}.nc
    if [[ ${model} == 'ICE' ]]; then
        file_out=${DTG_TEXT}.cice_model.res.nc
    else
        file_out=${DTG_TEXT}.${f}.nc
    fi
    WGET_AWS ${aws_path_replay}/${file_in} ${file_out} 
done
if [[ ${model} == 'ICE' ]]; then
    ${SCRIPT_DIR}/CICE_ic_edit.py -f ${file_out}
    echo "MOVING ${file_out%.nc}_new.nc" "${file_out}"
    mv "${file_out%.nc}_new.nc" "${file_out}"
fi

echo "SUCCESSFULLY Downloaded ${model}"
exit 0


