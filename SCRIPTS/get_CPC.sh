#!/bin/bash
set -ux
export dtg=${1}
model=${2} # ATM, OCN, ICE
export IC_SRC=${3:-CPC}
DEBUG=${DEBUG:-F}
source ${SCRIPT_DIR}/functions.sh
source ${SCRIPT_DIR}/defaults.sh

# glore path based on date
if [[ "${dtg}" -le 2004123100 ]]; then
    glore_c5path=/gpfs/f5/cfsrl/scratch/JieShun.Zhu/ng-godas/EXPrt.ice/BIAScorr.noARCTIC_CPCdeaggrT/rst
elif [[ "${dtg}" -ge 2005010100 ]] && [[ "${dtg}" -le 2021123100 ]]; then
    glore_c5path=/gpfs/f5/cfsrl/scratch/JieShun.Zhu/ng-godas/EXPrt.ice/BIAScorr.noARCTIC_CPCdeaggrTfnmoc2/rst
elif [[ "${dtg}" -ge 2022010100 ]] && [[ "${dtg}" -le 2024080100 ]]; then
    glore_c5path=/gpfs/f5/cpchso/scratch/JieShun.Zhu/ng-godas_1/EXPrt.ice/CALbump/rst
elif [[ "${dtg}" -ge 2026010100 ]]; then
    glore_c5path=/gpfs/f5/cpcice/scratch/JieShun.Zhu/ng-godas/EXPrt.ice/CALbump2/ana
fi

# core data path
if [[ "${dtg}" -le 2023123100 ]]; then
    core_path=/Permanent/NCEPDEV/cpc-om/Wesley.Ebisuzaki/core/nemsio/ens_00Z
elif [[ "${dtg}" -ge 2024010200 ]] && [[ "${dtg}" -le 2025050100 ]]; then
    core_path=/NCEPDEV/cpc-om/1year/Leigh.Zhang/core/flux
elif [[ "${dtg}" -ge 202506010100 ]]; then
    #The data, from 202507 to 202511, have been archived on HPSS at 
    core_path=/NCEPDEV/cpc-om/Permanent/Leigh.Zhang/core/flux
fi

# get file names
if [[ ${model} == OCN ]]; then
    file_in=${glore_c5path}/${dtg:0:8}12/ctrl/MOM.res.nc 
    file_out=${dir_restart_ocean}/${DTG_TEXT}.MOM.res.nc
    mkdir -p ${dir_restart_ocean}
elif [[ ${model} == ICE ]]; then
    file_in=${glore_c5path}/${dtg:0:8}12/ctrl/iced.${dtg:0:4}-${dtg:4:2}-${dtg:6:2}-43200.nc 
    file_out=${dir_restart_ice}/${DTG_TEXT}.cice_model.res.nc
    mkdir -p ${dir_restart_ice}
elif [[ ${model} == ATM ]]; then
    htar_file=${core_path}/ens_nem_${dtg:0:6}.tar
fi

# download
if [[ ${model} == OCN ]] || [[ ${model} == ICE ]]; then
    UUID_C5="57da978b-4434-4ac0-ae5e-3de5161305e1"
    if [[ ${machine} == "ursa" ]]; then
        UUID_DES="8bbc822a-a625-4058-b68e-65edacd73828"
        module load globus-cli
    else
        echo "Set Up globus SRC ID"
        echo "  globus endpoint search "NOAA RDHPCS ${machine}""
        exit 1
    fi

    TRANSFER_OUT=$(globus transfer "${UUID_C5}:${file_in}" "${UUID_DES}:${file_out}")
    TASK_ID=$(echo "$TRANSFER_OUT" | tail -n 1 | awk '{print $3}')
    
    # Safeguard: Check if TASK_ID is empty
    if [[ -z "${TASK_ID}" ]]; then
        echo "ERROR: Globus transfer failed to start. Output was:"
        echo "$TRANSFER_OUT"
        exit 1
    fi
    #TASK_ID=$( globus transfer ${UUID_C5}:${file_in} ${UUID_DES}:${file_out} | tail -n 1 | awk '{print $3}' )

    globus task wait -H --timeout 60  ${TASK_ID}
# ATM
else 
    mkdir -p ${dir_restart_atmos}
    cd ${dir_restart_atmos}
    htar -xvf ${htar_file}
    [[ ${?} > 0 ]] && echo "FATAL: htar failed" && exit 1
    # remove all files not associated with dtg
    find . -type f ! -name "*${dtg}*mem012*" -delete
    files=$(ls *gz)
    for f in ${files}; do
        gunzip ${f}
    done
fi

exit 0
