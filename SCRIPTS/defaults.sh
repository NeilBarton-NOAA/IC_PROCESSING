set -u
############
# defaults
SCRIPT_DIR=${SCRIPT_DIR:-$PWD}
IC_SRC=${IC_SRC:-"SCOUT"} #CPC, GFS
ATMRES=${ATMES:-"C192"}
OCNRES=${OCN:-"mx025"}
run=${run:-sfs}
GLOBUS=F
dtg_minus6=$(date -u -d"${dtg:0:4}-${dtg:4:2}-${dtg:6:2} ${dtg:8:2}:00:00 6 hours ago" +%Y%m%d%H)
dtg_minus3=$(date -u -d"${dtg:0:4}-${dtg:4:2}-${dtg:6:2} ${dtg:8:2}:00:00 3 hours ago" +%Y%m%d%H)
dtg_plus3=$(date -u -d"${dtg:0:4}-${dtg:4:2}-${dtg:6:2} ${dtg:8:2}:00:00 3 hours" +%Y%m%d%H)
NENS=10

############
# scrath dir based in machine
machine=$(uname -n)
export chgres_compiler=intel
case ${machine} in 
    gaea* | dtn* | c6*)
        SCRATCH_DIR=/gpfs/f6/sfs-emc/scratch/${USER}
        m_target=gaeac6
        ;;
    hercules*)
        SCRATCH_DIR=/work/noaa/marine/${USER}
        m_target=hercules
        ;;
    ufe* | u*c*)
        SCRATCH_DIR=/scratch4/NCEPDEV/stmp/${USER}
        m_target=ursa
        export chgres_compiler=intelllvm
        ;;
    *)
        echo "FATAL: machine unknown: ${machine}" && exit 1
        ;;
esac
CODE_DIR=${SCRATCH_DIR}/CODE/IC_PROCESSING
if [[ ${IC_SRC} == "SCOUT" ]]; then
    SCOUT_VERSION=${SCOUT_VERSION:-'SFv1.01'}
    export IC_DIR=${SCRATCH_DIR}/ICs/${SCOUT_VERSION}/${ATMRES}${OCNRES} 
elif [[ ${IC_SRC} == "GFS" ]]; then
    export IC_DIR=${SCRATCH_DIR}/ICs/GFS/${ATMRES}${OCNRES} 
else
    echo "FATAL: IC_SRC not set up", ${IC_SRC} && exit 1
fi
export IC_DIR=${ICDIR:-$IC_DIR} && mkdir -p ${IC_DIR}

###########
# download directories for each component
export dir_restart_atmos=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/atmos/restart
export dir_restart_ocean=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/ocean/restart
export dir_restart_ice=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/ice/restart
export dir_restart_wave=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/wave/restart
export dir_restart_med=${IC_DIR}/${run}.${dtg_minus6:0:8}/${dtg_minus6:8:2}/mem000/model/med/restart
export dir_input_atmos=${IC_DIR}/${run}.${dtg:0:8}/${dtg:8:2}/mem000/model/atmos/input
export dir_inc_atmos=${IC_DIR}/${run}.${dtg:0:8}/${dtg:8:2}/mem000/analysis/atmos
export dir_inc_ocean=${IC_DIR}/${run}.${dtg:0:8}/${dtg:8:2}/mem000/analysis/ocean

############
# files
export restart_tile_files_atmos='ca_data fv_core.res fv_srf_wnd.res fv_tracer.res phy_data sfc_data'
export restart_nontile_files_atmos='ca_data fv_core.res'

############
# GFS Retro Run
#   Aug 30 2022 to Oct 10 2022
#   Mar  1 2024 to Nov 30 2025
#   Nov 20 2025 to Feb 28 2026 (near real time run)
#   ICs may be on Mondays
export hpss_path=/5year/NCEPDEV/emc-global/emc.glopara/*/GFSv17/retrov17*

############
# CPC ICs
hpss_core_path=/Permanent/NCEPDEV/cpc-om/Wesley.Ebisuzaki/core/nemsio/ens_00Z
glore_path=onc5 #(uce scp command)

############
# SFS Scout Run
# https://noaa-reanalyses-pds.s3.amazonaws.com/index.html
export aws_path="https://noaa-reanalyses-pds.s3.amazonaws.com/analyses/scout_runs/3dvar_coupledreanl_scoutrun_v1.01"
export DTG_TEXT_SRC=${dtg_plus3:0:8}.${dtg_plus3:8:10}0000 
export DTG_TEXT_DES=${dtg:0:8}.${dtg:8:10}0000 
export DTG_TEXT=${DTG_TEXT_DES}
########################
# GEFSv13 Replay
# https://noaa-ufs-gefsv13replay-pds.s3.amazonaws.com/index.html
# aws_path="https://noaa-ufs-gefsv13replay-pds.s3.amazonaws.com/${dtg:0:4}/${dtg:4:2}/${dtg:0:8}06"
# aws_C192sfc="https://noaa-oar-sfsdev-pds.s3.amazonaws.com/input/c192/hr4_land/${dtg}"

########################
# for chgres
export APRUN="srun"

