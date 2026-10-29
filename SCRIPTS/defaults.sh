set -u
############
# defaults
SCRIPT_DIR=${SCRIPT_DIR:-$PWD}
IC_SRC=${IC_SRC:-"GFS"} # REPLAY, GFS, SCOUT,GFS_DA
ATMRES=${ATMRES:-"C192"}
OCNRES=${OCNRES:-"mx025"}
run=${run:-sfs}
NENS=${NENS:-10}
source ${SCRIPT_DIR}/../MACHINE/config.sh

# location of code and downloads
if [[ ${IC_SRC} == "SCOUT" ]]; then
    SCOUT_VERSION=${SCOUT_VERSION:-'SFv1.01'}
    export IC_DIR=${WORK_DIR}/ICs/${SCOUT_VERSION}/${ATMRES}${OCNRES} 
else
    export IC_DIR=${WORK_DIR}/ICs/${IC_SRC}/${ATMRES}${OCNRES} 
fi
export IC_DIR=${ICDIR:-$IC_DIR} && mkdir -p ${IC_DIR}

# dtg specific
dtg_minus6=$(date -u -d"${dtg:0:4}-${dtg:4:2}-${dtg:6:2} ${dtg:8:2}:00:00 6 hours ago" +%Y%m%d%H)
dtg_minus3=$(date -u -d"${dtg:0:4}-${dtg:4:2}-${dtg:6:2} ${dtg:8:2}:00:00 3 hours ago" +%Y%m%d%H)
dtg_plus3=$(date -u -d"${dtg:0:4}-${dtg:4:2}-${dtg:6:2} ${dtg:8:2}:00:00 3 hours" +%Y%m%d%H)

# file prefix
export DTG_TEXT=${dtg:0:8}.${dtg:8:10}0000 
export DTG_TEXT_PLUS3=${dtg_plus3:0:8}.${dtg_plus3:8:10}0000 
export DTG_TEXT_MINUS3=${dtg_minus3:0:8}.${dtg_minus3:8:10}0000 

# download directories for each component
export dir_restart_atmos=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/atmos/restart
export dir_restart_ocean=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/ocean/restart
export dir_restart_ice=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/ice/restart
export dir_restart_wave=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/wave/restart
export dir_restart_med=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/med/restart
export dir_input_atmos=${IC_DIR}/${run}.${dtg:0:8}/${dtg:8:2}/mem000/model/atmos/input
export dir_inc_atmos=${IC_DIR}/${run}.${dtg:0:8}/${dtg:8:2}/mem000/analysis/atmos
export dir_inc_ocean=${IC_DIR}/${run}.${dtg:0:8}/${dtg:8:2}/mem000/analysis/ocean
export dir_inc_ice=${IC_DIR}/${run}.${dtg:0:8}/${dtg:8:2}/mem000/analysis/ice

# restart files
export restart_tile_files_atmos='ca_data fv_core.res fv_srf_wnd.res fv_tracer.res phy_data sfc_data'
export restart_nontile_files_atmos='ca_data fv_core.res'
export restart_files_ocean='MOM.res MOM.res_1 MOM.res_2 MOM.res_3'

# for chgres
m_target=${machine}

