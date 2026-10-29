#!/bin/bash
set -xu
dtg=${1}
IC_SRC=${2}
SCRIPT_DIR=${SCRIPT_DIR:-$(dirname "$0")}
source ${SCRIPT_DIR}/defaults.sh

LINK_MEMBERS () {
dir=${1}
NENS=${2}
mem000_dir="../../../../mem${dir##*mem}"
for i in $(seq 1 ${NENS}); do
    mem=$(printf "%03d" ${i})   
    dir_mem=${dir/mem000/mem${mem}} 
    mkdir -p ${dir_mem} && cd ${dir_mem}
    echo $dir_mem
    files=$( ls ${mem000_dir}/* )
    for f in ${files}; do
        ln -sf ${mem000_dir}/$(basename ${f}) .
    done
done
}

#For linking ice analysis for GFS_DA ICs
LINK_MEMBERS2 () {
dir=${1}
NENS=${2}
    mem000_dir="../../../mem${dir##*mem}"
for i in $(seq 1 ${NENS}); do
    mem=$(printf "%03d" ${i})
    dir_mem=${dir/mem000/mem${mem}}
    mkdir -p ${dir_mem} && cd ${dir_mem}
    echo $dir_mem
    files=$( ls ${mem000_dir}/* )
    for f in ${files}; do
        ln -sf ${mem000_dir}/$(basename ${f}) .
    done
done
}

[[ -d ${dir_restart_atmos} ]] && LINK_MEMBERS ${dir_restart_atmos} ${NENS}
[[ -d ${dir_input_atmos} ]] && LINK_MEMBERS ${dir_input_atmos} ${NENS}
[[ -d ${dir_restart_ocean} ]] && LINK_MEMBERS ${dir_restart_ocean} ${NENS}
[[ -d ${dir_restart_ice} ]] && LINK_MEMBERS ${dir_restart_ice} ${NENS}
[[ -d ${dir_restart_wave} ]] && LINK_MEMBERS ${dir_restart_wave} ${NENS}
[[ -d ${dir_restart_med} ]] && LINK_MEMBERS ${dir_restart_med} ${NENS}
[[ -d ${dir_inc_ice} ]] && LINK_MEMBERS2 ${dir_inc_ice} ${NENS}

echo "Linked ICs to member directories"
