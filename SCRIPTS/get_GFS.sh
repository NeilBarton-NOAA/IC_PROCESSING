#!/bin/bash
set -u
export dtg=${1}
file=${2} # gdasocean_restart, gdas_restartb, gdasice_restart, enkfgdas_restartb_grp1, enkfgdas_restarta_grp1
DEBUG=${DEBUG:-F}
#SCRIPT_DIR=${SCRIPT_DIR:-$(dirname "$0")}
export IC_SRC=GFS
source ${SCRIPT_DIR}/functions.sh
source ${SCRIPT_DIR}/defaults.sh
dir=${IC_DIR}
PREFIX="SFS"
f_extracted=${dir}/${dtg}_${file}_htar.log

echo "DOWNLOADING GFS RESTARTS to ${dir}"
mkdir -p ${dir} && cd ${dir}
if [[ ${DOWNLOAD:-T} == T ]]; then
    GFS_RESTART_DTG ${dtg} ${hpss_path} 
    echo "Target: ${dtg}" 
    echo "Closest match: ${dtg_closest}" 
    echo "Days Apart: ${day_diff}" 
    dtg_closest_minus6=$(date -d"${dtg_closest:0:8} ${dtg_closest:8:2} 6 hours ago" +%Y%m%d%H)
    echo $file

    if [[ "${file}" == "gdasocean_restart" ]] || [[ "${file}" == *"restartb"* ]]; then
        hpss_file="${hpss_file//${dtg_closest}/${dtg_closest_minus6}}"    
    fi
    hpss_file="${hpss_file//enkfgdas_restarta_grp1.tar/${file}.tar}"    
    echo "Downloading: ${hpss_file}"
    [[ -f ${f_extracted} ]] && rm ${f_extracted}
    htar -xvf ${hpss_file} > ${f_extracted} 2>&1
    if (( ${?} > 0 )); then
        echo 'FATAL in htar, file also not at'
        echo '  hpss_file:', ${hpss_file}
        exit 1
    fi
else
    dtg_closest=2025062900
    dtg_closest_minus6=$(date -d"${dtg_closest:0:8} ${dtg_closest:8:2} 6 hours ago" +%Y%m%d%H)
fi
############
# rename to dtg_minus
if [[ ${file} == 'gdasocean_restart' ]]; then
    files=$( grep MOM ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | sort  )
    for f in ${files}; do
        MV ${f} ${PREFIX}
    done        
elif [[ ${file} == 'gdasocean_analysis' ]]; then
    f=$(grep mom6_incre ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
    MV ${f} ${PREFIX}
    f=$(grep cice ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
    MV ${f} ${PREFIX}
elif [[ ${file} == 'gdas_restarta' ]]; then
    sfc_files="sfcanl_data increment.sfc"
    for f_res in ${sfc_files}; do
        for t in {1..6}; do
            f=$(grep ${f_res}.tile${t}.nc ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
            MV ${f} ${PREFIX}
        done
    done
    files=$(grep gdas.t00z.increment.atm ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1)
    for f in ${files}; do
        MV ${f} ${PREFIX}
    done
elif [[ ${file} == 'gdas_restartb' ]]; then
    for t in {1..6}; do
        files='ca_data fv_core.res fv_srf_wnd.res fv_tracer.res phy_data'
        for f_res in ${restart_tile_files_atmos}; do
            f=$(grep ${dtg_closest_minus6:0:8}.210000.${f_res}.tile${t}.nc ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
            MV ${f} ${PREFIX}
        done        
    done
    #fv_core.res
    f=$(grep ${dtg_closest_minus6:0:8}.210000.fv_core.res.nc ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
    MV ${f} ${PREFIX}
elif [[ ${file} == 'enkfgdas_restarta_grp1' ]]; then
    members=$( grep 210000.analysis.cice_model.res ${f_extracted} | cut -d' ' -f3 | cut -d'/' -f3 )
    for mem in ${members}; do
        # cice restart
        f=$(grep cice_model ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
        MV ${f} ${PREFIX}
        # mom6 increment
        f=$(grep mom6_incre ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
        MV ${f} ${PREFIX}
        # atmos non tiles
        files=$( grep atmos ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 )
        for f in ${files}; do
            MV ${f} ${PREFIX}
        done
        exit 1
    done
elif [[ ${file} == 'enkfgdas_restartb_grp1' ]]; then
    members=$( grep 210000.fv_core.res.tile1 ${f_extracted} | cut -d' ' -f3 | cut -d'/' -f3 )
    for mem in ${members}; do
        # atmos files
        for f_res in ${restart_tile_files_atmos}; do
            for t in {1..6}; do
                f=$(grep ${dtg_closest_minus6:0:8}.210000.${f_res}.tile${t}.nc ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
                MV ${f} ${PREFIX}
            done
        done
        #fv_core.res
        f=$(grep ${dtg_closest_minus6:0:8}.210000.fv_core.res.nc ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
        MV ${f} ${PREFIX}
        # mom6 files
        f=$(grep ${dtg_closest_minus6:0:8}.210000.MOM.res.nc ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
        MV ${f} ${PREFIX}
        for t in {1..3}; do
            f=$(grep ${dtg_closest_minus6:0:8}.210000.MOM.res_${t}.nc ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
            MV ${f} ${PREFIX}
        done 
    done
else
    echo "FATAL Script is not set up for ${file}"
    exit 1
fi

#if [[ ${DEBUG:-F} == F ]]; then
#    rm -r ${f%%/*}
#else
#    echo "NPB Check" && exit 1
#fi
echo "Downloaded ${file}"
exit 0


