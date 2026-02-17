#!/bin/sh
set -u
WGET_AWS () {
file_in=${1}
file_out=${2}
if [[ ! -s ${file_out} ]]; then
    rm ${file_out}
fi
if [[ ! -f ${file_out} ]]; then
    wget ${file_in} -O ${file_out} 
    echo "wget ${file_in} -O ${file_out}"
else
    echo "${file_out} already exist"
fi
}

GLOBUS_AWS () {
file_in=${1}
file_out=${2}
if [[ ! -f ${file_out} ]]; then
    ID=$( globus transfer ${UUID_AWS_S3_PUBLIC}\://${file_in} \
        ${UUID_HERA_DTN}\:${file_out} | \
        tail -n 1 | \
        awk '{print $3}' )
    (( ${#ID} != 36 )) && ID=9999
    echo ${ID}
fi
}

FIND_EMPTY_FILES () {
dir_in=${1}
n_empty=$( find ${dir_in} -type f -size -17k | grep -v ca_data | grep -v coupler.res | wc -l )
if (( ${n_empty} >> 0 )); then
    echo "Failed: empty files found"
    files=$( find ${dir_in} -type f -size -17k | grep -v ca_data | grep -v coupler.res )
    for f in ${files}; do
        echo "  removing:  "${f}
        rm ${f}
    done
    #find ${dir_in} -type f -size 0 -delete
    exit 1
fi
}

FIND_INDEX() {
local target="$1"
shift
local arr=("$@")
for i in "${!arr[@]}"; do
    if [[ "${arr[$i]}" == "$target" ]]; then
        echo "$i"
        return 0
    fi
done
return 1
}

RENAME_SFS() {
    local f=${1}
    local nf=${f}
    nf="${nf//enkfgdas/sfs}"
    nf="${nf//gdas/sfs}"
    nf="${nf//${dtg_closest_minus6:0:8}/${dtg_minus6:0:8}}"  
    nf="${nf//${dtg_closest:0:8}/${dtg:0:8}}"  
    if [[ ! ${nf} == *"mem"* ]]; then
        nf=$(echo "${nf}" | sed 's|\(/[0-9][0-9]/\)|\1mem000/|')
    fi
    echo ${nf}
}

RENAME_GFS() {
    local f=${1}
    local nf="${f//${dtg_closest:0:8}/${dtg_minus6:0:8}}"    
    nf="${nf//${dtg_closest_plus6:0:8}/${dtg:0:8}}"  
    echo ${nf}
}

MV() {
    f1=${1}
    local PREFIX=${2}
    if [[ ${PREFIX} == "SFS" ]]; then
        f2=$( RENAME_SFS ${f1} )
    else
        f2=$( RENAME_GFS ${f1} )
    fi
    mkdir -p $( dirname ${f2} )
    if [[ ${MV_DATA:-"T"} == "T" ]]; then
        echo "mv ${f1} ${f2} ${PREFIX}"
        mv ${f1} ${f2}
    else
        echo "ln -sf ${f1} ${f2} ${PREFIX}"
        f1=$(realpath --relative-to=$(dirname "$f2") "$f1")
        ln -sf ${f1} ${f2}
    fi
}

GFS_RESTART_DTG() {
    dtg=${1}
    hpss_path=${2}
    files=() && dtgs=()
    all_files=$( hsi find ${hpss_path} -name "enkfgdas_restarta_grp1.tar" 2>&1 | grep NCEPDEV )
    for f in ${all_files}; do
        files+=(${f})
        dtgs+=($(echo "${f}" | cut -d'/' -f9))
    done
    sec_target=$(date -d "${dtg:0:8}" +%s)
    dtg_closest="" && min_diff=-1
    for d in ${dtgs[@]}; do
        # Convert current date in loop to seconds
        sec_current=$(date -d "${d:0:8}" +%s)
        # Calculate absolute difference
        diff=$(( sec_current - sec_target ))
        abs_diff=${diff#-} # This removes the negative sign if it exists
        # Check if this is the smallest difference found so far
        if [[ ${min_diff} -eq -1 ]] || [[ ${abs_diff} -lt ${min_diff} ]]; then
            min_diff=${abs_diff}
            dtg_closest=${d}
        fi
    done
    day_diff=$(( min_diff / 86400 ))
    if [[ ${day_diff} > 7 ]]; then
        echo "FATAL Closest Match is Too Far Away from Target:" ${day_diff} && exit 1
    fi 
    for i in "${!dtgs[@]}"; do
        if [[ "${dtgs[${i}]}" == "${dtg_closest}" ]]; then
            hpss_file=${files[${i}]}
            return 0
        fi
    done
}

