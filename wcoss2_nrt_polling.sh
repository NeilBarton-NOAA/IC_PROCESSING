#!/bin/bash
set -u
dtg=${1:-2026110100}
dir_nrt="/lfs/h2/emc/gfstemp/emc.global/comroot/retrov17_01_realtime"
HOME_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
SCRIPT_DIR=${HOME_DIR}/SCRIPTS
export COMPUTE_ACCOUNT=GFS-DEV
source ${SCRIPT_DIR}/defaults.sh
dir_sfs=${WORK_DIR}/ICs/GFS/C192mx025
main(){
local dtg=$1
echo ""
current_date=$(date)
echo "CURRENT date ${current_date}"
echo    'dir of realtime GFS', ${dir_nrt}
echo    'dir of sfs ice',      ${dir_sfs}
n_files=$( find ${dir_nrt}/*.${dtg:0:8}/${dtg:8:10} -name *mom6_increment*.nc 2>/dev/null | wc -l )
if [[ ${n_files} == 0 ]]; then
    echo "FILES FOR ${dtg} not found"
    exit 0
fi

COUNT_FILES
if (( ${NUM_FILES} == 1605 )); then
    echo "FILES ALREADY PROCESSED FOR" ${dtg}
    exit 0
fi

FIND_AND_COPY "${restart_tile_files_atmos} analysis.cice_model.res ${restart_files_ocean}"
FIND_AND_COPY "increment.sfc increment.atm mom6_increment ensmean_increment.sfc recentered_increment.atm"

${HOME_DIR}/RUN_CHGRES.sh ${dtg} GFS
${HOME_DIR}/RUN_REGRID.sh ${dtg} GFS
while [[ $(qstat -u "${USER}" | grep -E -c "CHGRES|REGRID") -gt 0 ]]; do
    echo "WAITING UNTIL CHGRES and REGRID Jobs are Finished"
    sleep 60  # Wait 60 seconds before checking again
done

# htar into groups
HTAR_MEMBERS 000 010
HTAR_MEMBERS 011 020
HTAR_MEMBERS 021 030
}

COUNT_FILES(){
dir1=${dir_sfs}/sfs.${dtg:0:8}
dir2=${dir_sfs}/sfs.${dtg_minus6:0:8}
NUM_FILES=$(find ${dir1} ${dir2} -type f 2>/dev/null | wc -l)
echo "NUMBER of FILES in ${dir1} and ${dir2}: $NUM_FILES"
}

FIND_AND_COPY(){
local restart_files=${1}
for rf in ${restart_files}; do
    if [[ ${rf} == *cice* ]]; then
        DTG_TEXT=${DTG_TEXT_MINUS3}
        dtg_dir=${dtg}
    elif [[ ${rf} == *increment* ]]; then
        DTG_TEXT=""
        dtg_dir=${dtg}
    else
        DTG_TEXT=${DTG_TEXT_MINUS3}
        dtg_dir=${dtg_minus6}
    fi
    RUN="*"
    if [[ ${rf} == increment.atm ]]; then
        RUN="gdas"
    fi
    echo ${rf}, ${dtg}, ${dtg_dir}, ${DTG_TEXT}
    files=$( find ${dir_nrt}/${RUN}.${dtg_dir:0:8}/${dtg_dir:8:10} -name *${DTG_TEXT}.${rf}*.nc | sort )
    for f in ${files}; do
        RENAME_COPY ${f}
    done
done

}

RENAME_COPY() {
    local f=${1}
    local nf=${f}
    nf="${nf//enkfgdas/sfs}"
    nf="${nf//gdas/sfs}"
    nf="${nf//${dir_nrt}/${dir_sfs}}"
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
    ############
    # COPY
    mem=$( echo ${f2} | grep -oP '(?<=mem)\d+' )
    [[ -f ${f2} ]] && echo "FILE ALREADY COPIED: ${f2}"
    if (( 10#${mem} < 31 )) && [[ ${f2} != *gfs* ]] && [[ ! -f ${f2} ]]; then
    #if (( 10#${mem} < 3 )) && [[ ${f2} != *gfs* ]] && [[ ! -f ${f2} ]]; then
        mkdir -p $(dirname ${f2})
        cp ${f} ${f2} 
        echo "COPIED: ${f2}"
    fi
}

HTAR_MEMBERS() {
    mem_l=${1}
    mem_h=${2}
    f=/NCEPDEV/emc-marine/5year/Neil.Barton/SFS_GFS_ICS/${dtg}_C192mx025_mem${mem_l}_to_mem${mem_h}.tar
    ds=""
    for mem in $(seq ${mem_l} ${mem_h}); do
        m=$(printf "%03d" ${mem})
        ds="${ds} sfs.${dtg:0:8}/${dtg:8:2}/mem${m}/* sfs.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem${m}/*"
    done
    JOB_NAME=HTAR.SFSICS.${dtg}.${mem_l}to${mem_h}
    WALLTIME="03:00:00"
    source ${HOME_DIR}/MACHINE/config.sh
    echo ${JOB_NAME}
    echo -e "${SUBMIT_HPSS}\ncd ${dir_sfs}\nhtar -cvf ${f} ${ds}"> submit_HTAR.sh
    qsub submit_HTAR.sh
    rm submit_HTAR.sh
}

main ${dtg}
[[ "${USER}" == "neil.barton" ]] && rocotoDEAD -r
