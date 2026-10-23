#!/bin/bash
set -u
set -x
export dtg=${1}
file=${2} # gdasocean_restart, gdas_restartb, gdasice_restart, enkfgdas_restartb_grp1, enkfgdas_restarta_grp1
DEBUG=${DEBUG:-F}
export IC_SRC=${3:-GFS}
source ${SCRIPT_DIR}/functions.sh
source ${SCRIPT_DIR}/defaults.sh
dir=${IC_DIR}
f_extracted=${dir}/${dtg}_${file}_htar.log

############
# GFS Retro Run
#   Aug 30 2022 to Oct 10 2022
#   Mar  1 2024 to Nov 30 2025
#   Nov 20 2025 to Feb 28 2026 (near real time run)
#   ICs may be on Mondays
if [[ ${IC_SRC} == "GFS" ]]; then
    export hpss_path=/5year/NCEPDEV/emc-global/emc.glopara/*/GFSv17/retrov17*
elif [[ ${IC_SRC} == "GFS_DA" ]]; then
   if [[ "${dtg}" -ge 2025091512 ]] && [[ "${dtg}" -le 2026042612 ]]; then
       export hpss_path="/NCEPDEV/emc-global/1year/john.steffen/WCOSS2/scratch/cp4.03-parallel-hybrid" # 2025091512 -
   elif [[ "${dtg}" -ge 2024010312 ]] && [[ "${dtg}" -le 2024080612 ]]; then
       export hpss_path="/NCEPDEV/emc-global/5year/emc.da/WCOSS2/scratch/gfsv17-marine3dvar-reanalysis-stream2" # 2024010312 - 2024080612
   elif [[ "${dtg}" -ge 2021070212 ]] && [[ "${dtg}" -le 2022042618 ]]; then
       export hpss_path="/NCEPDEV/emc-global/5year/emc.da/WCOSS2/scratch/gfsv17-marine3dvar-reanalysis" # 2021070212 - 2022042618
   fi
fi

# download 
echo "GFS RESTARTS to ${dir}"
echo "  ${f_extracted}"
mkdir -p ${dir} && cd ${dir}
GFS_RESTART_DTG ${dtg} ${hpss_path} gdasocean_restart.tar
echo "Target: ${dtg}" 
echo "Closest match: ${dtg_closest}" 
echo "Days Apart: ${day_diff}" 
hpss_file=${hpss_file//gdasocean_restart/${file}}
dtg_closest_minus6=$(date -d"${dtg_closest:0:8} ${dtg_closest:8:2} 6 hours ago" +%Y%m%d%H)
if [[ ${file} =~ ^(gdas_restarta|enkfgdas_restarta_grp1|gdasocean_analysis|enkfgdas)$ ]] ; then
    hpss_file="${hpss_file//${dtg_closest_minus6}/${dtg_closest}}"
fi
if [[ ${DOWNLOAD:-T} == T ]]; then   
    echo "Downloading: ${hpss_file}"
    [[ -f ${f_extracted} ]] && rm ${f_extracted}
    htar -xvf ${hpss_file} > ${f_extracted} 2>&1
    if (( ${?} > 0 )); then
        echo 'FATAL in htar, file also not at'
        echo '  hpss_file:', ${hpss_file}
        exit 1
    fi
fi

# rename to dtg_minus
if [[ "${file}" == "gdasocean_restart" ]]; then
    files=$( grep MOM ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | sort  )
    [[ "${#files}" == 0 ]] && exit 1
    for f in ${files}; do
        MV ${f}
        [[ $? > 0 ]] && exit 1
    done        
elif [[ "${file}" == "gdasocean_analysis" ]]; then
    f=$(grep mom6_incre ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
    MV ${f}
    [[ $? > 0 ]] && exit 1
    f=$(grep cice ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
    MV ${f}
    [[ $? > 0 ]] && exit 1
elif [[ "${file}" == "gdas_restarta" ]]; then
    if [[ ${IC_SRC} == "GFS" ]]; then
        sfc_files="sfcanl_data increment.sfc"
    elif [[ ${IC_SRC} == "GFS_DA" ]]; then
        sfc_files="sfcanl_data"
    fi
    for f_res in ${sfc_files}; do
        for t in {1..6}; do
            f=$(grep ${f_res}.tile${t}.nc ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
            MV ${f}
            [[ $? > 0 ]] && exit 1
        done
    done
    files=$(grep gdas.t00z.increment.atm ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1)
    [[ "${#files}" == 0 ]] && exit 1
    for f in ${files}; do
        MV ${f}
        [[ $? > 0 ]] && exit 1
    done
elif [[ "${file}" == "gdas_restartb" ]]; then
    for t in {1..6}; do
        files='ca_data fv_core.res fv_srf_wnd.res fv_tracer.res phy_data'
        for f_res in ${restart_tile_files_atmos}; do
            f=$(grep "${dtg_closest_minus6:0:8}.210000.${f_res}.tile${t}.nc" ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
            MV ${f}
            [[ $? > 0 ]] && exit 1
        done        
    done
    f=$(grep "${dtg_closest_minus6:0:8}.210000.fv_core.res.nc" ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
    MV ${f}
    [[ $? > 0 ]] && exit 1
elif [[ "${file}" == "enkfgdas" ]]; then
    files=$( grep ensmean_increment.sfc.i ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 )
    [[ "${#files}" == 0 ]] && exit 1
    for f in ${files}; do
        MV ${f}
        [[ $? > 0 ]] && exit 1
    done
elif [[ "${file}" == "enkfgdas_restarta_grp"* ]]; then
    members=$( grep 210000.analysis.cice_model.res ${f_extracted} | cut -d' ' -f3 | cut -d'/' -f3 )
    for mem in ${members}; do
        # cice restart
        f=$(grep cice_model ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
        MV ${f}
        [[ $? > 0 ]] && exit 1
        # increments mom6, atm, and sfc
        files=$(grep increment ${f_extracted} | grep -v tile | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 )
        for f in ${files}; do
            MV ${f}
            [[ $? > 0 ]] && exit 1
        done
    done
elif [[ "${file}" == "enkfgdas_restartb_grp"* ]]; then
    members=$( grep 210000.fv_core.res.tile1 ${f_extracted} | cut -d' ' -f3 | cut -d'/' -f3 )
    for mem in ${members}; do
        if [[ ${IC_SRC} != "GFS_DA" ]]; then
            # atmos files
            for f_res in ${restart_tile_files_atmos}; do
                for t in {1..6}; do
                    f=$(grep ${dtg_closest_minus6:0:8}.210000.${f_res}.tile${t}.nc ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
                    MV ${f}
                    [[ $? > 0 ]] && exit 1
                done
            done
            #fv_core.res
            f=$(grep ${dtg_closest_minus6:0:8}.210000.fv_core.res.nc ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
            MV ${f}
            [[ $? > 0 ]] && exit 1
        fi
        # mom6 files
        f=$(grep ${dtg_closest_minus6:0:8}.210000.MOM.res.nc ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
        MV ${f}
        [[ $? > 0 ]] && exit 1
        for t in {1..3}; do
            f=$(grep ${dtg_closest_minus6:0:8}.210000.MOM.res_${t}.nc ${f_extracted} | grep ${mem} | cut -d' ' -f3 | cut -d',' -f1 | head -n 1)
            MV ${f}
            [[ $? > 0 ]] && exit 1
        done 
    done
else
    echo "FATAL Script is not set up for ${file}"
    exit 1
fi

if [[ ${MV_DATA:-"T"} == "T" ]]; then
    files=$( grep 'HTAR: x' ${f_extracted} | cut -d' ' -f3 | cut -d',' -f1 )
    for f in ${files}; do
        [[ -f ${f} ]] && rm ${f}
    done
fi

echo "SUCCESSFULLY Downloaded ${file}"
exit 0
