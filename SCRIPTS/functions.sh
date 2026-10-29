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
    if [[ "${nf}" == *"/18/"* ]]; then
        nf="${nf//${dtg_closest_minus6:0:8}/${dtg_minus6:0:8}}"
    else
        nf="${nf//${dtg_closest:0:8}/${dtg:0:8}}"
    fi
    if [[ "${nf}" == *"cice"* ]]; then
        old="${dtg_closest_minus6:0:8}.210000"
        if [[ ${IC_SRC} == "GFS_DA" ]]; then
            new="${dtg:0:8}.000000"
        else
            new="${dtg_minus6:0:8}.210000"
        fi
        nf="${nf//${old}/${new}}"
    fi
        if [[ "${nf}" == *"MOM"* ]]; then
        old="${dtg_closest_minus6:0:8}.210000"
        if [[ ${IC_SRC} == "GFS_DA" ]]; then
            new="${dtg:0:8}.000000"
        else
            new="${dtg_minus6:0:8}.210000"
        fi
        nf="${nf//${old}/${new}}"
    fi
    if [[ ! ${nf} == *"mem"* ]]; then
        nf=$(echo "${nf}" | sed 's|\(/[0-9][0-9]/\)|\1mem000/|')
    fi
    if [[ ${nf} == *"increment.atm."* ]]; then
        # Extract the number
        num=$(echo "${nf}" | grep -oP '(?<=\.i)\d+(?=\.nc)')
        # Calculate new number
        new_num=$(echo "${num} - 3" | bc)
        # Replace in string
        nf=$(printf "${nf/.i${num}.nc/.i%03d.nc}" "${new_num}")
    fi
    f2=${nf}
}

MV() {
    f1=${1}
    RENAME_SFS ${f1}
    mkdir -p $( dirname ${f2} )
    if [[ ${MV_DATA:-"T"} == "T" ]]; then
        echo "mv ${f1} ${f2}"
        mv ${f1} ${f2}
    else
        echo "ln -sf ${f1} ${f2}"
        f1=$(realpath --relative-to=$(dirname "$f2") "$f1")
        ln -sf ${f1} ${f2}
    fi
}

GFS_RESTART_DTG() {
    dtg=${1}
    hpss_path=${2}
    local file=${3}
    local find_command=${4:-F}
    find_file="${SCRIPT_DIR}/../logs/${file}.log"
    
    if [[ ! -f ${find_file} ]] || [[ ${find_command} == T ]]; then
        hsi find "${hpss_path}" -name "${file}" 2>&1 | grep NCEPDEV | sort > "${find_file}"
    fi
    files=() && dtgs=()
    while IFS= read -r f; do
        files+=(${f})
        dtgs+=($(echo "${f}" | cut -d'/' -f9))
    done < ${find_file}
    sec_target=$(date -d "${dtg:0:4}-${dtg:4:2}-${dtg:6:2} ${dtg:8:2}:00:00 6 hours ago" +%s)
    dtg_closest="" && min_diff=-1
    for d in ${dtgs[@]}; do
        # Convert current date in loop to seconds
        sec_current=$(date -d "${d:0:4}-${d:4:2}-${d:6:2} ${d:8:2}:00:00" +%s)
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
            dtg_closest=$(date -d"${dtg_closest:0:8} ${dtg_closest:8:2} 6 hours" +%Y%m%d%H)
            return 0
        fi
    done
}

GET_AWS_TILES() {
files=${1}
}
