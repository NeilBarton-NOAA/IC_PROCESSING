#!/bin/bash
set -u
dtg=${1}
source ${SCRIPT_DIR}/defaults.sh
cd ${IC_DIR}

HPSS_DIR=/NCEPDEV/emc-marine/2year/${USER}/ICS/${IC_SRC}/${ATMRES}${OCNRES} 
folders=""
[[ -d ${IC_DIR}/${run}.${dtg:0:8} ]] && folders="${folders} ${run}.${dtg:0:8}"
[[ -d ${IC_DIR}/${run}.${dtg_minus6:0:8} ]] && folders="${folders} ${run}.${dtg_minus6:0:8}"

# Add logs if from GFS to see what day they are from
if [[ ${IC_SRC} == "GFS" ]]; then 
    files=$( ls ${dtg}*log )
    folders="${folders} ${files}"
fi

hsi mkdir -p ${HPSS_DIR}
file_name=${HPSS_DIR}/${dtg}.tar
htar -cvf ${file_name} ${folders}
if (( ${?} > 0 )); then
    echo 'ERROR in htar'
    exit 1
fi

if [[ ${IC_SRC} == "GFS" ]]; then 
    ds=$( grep 'HTAR: x' ${dtg}*log | cut -d' ' -f3 | cut -d'/' -f1 | sort | uniq )
    for d in ${ds}; do
     rm -r ${d}
    done
fi
echo "SUCCESFUL hpss_put_ICs.sh"
