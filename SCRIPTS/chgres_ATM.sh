#!/bin/sh
# https://noaa-emcufs-utils.readthedocs.io/en/latest/
set -u
dtg=${1}
SRC_ATMRES=${2:-"C192"}
SRC_OCNRES=${3:-"mx025"}
MEM=${4:-'000'}
source ${SCRIPT_DIR}/defaults.sh
DTG_TEXT=${dtg_minus3:0:8}.${dtg_minus3:8:10}0000 
compiler=${chgres_compiler}
dir=${dir_restart_atmos}
dir=${dir//mem000/mem${MEM}}
export DATA=${dir}
cd ${DATA}
echo ${DATA}
########################
# chgres_cube.sh options
export HOMEufs=${CODE_DIR}/UFS_UTILS
# needed for IC_SRC == "restart"
ORO_SRC=() && CORE_SRC=() && TRACER_SRC=() && SFC_SRC=() && ORO_DES=()
for i in {1..6}; do
    ORO_SRC+=("${SRC_ATMRES}.${SRC_OCNRES}_oro_data.tile${i}.nc")
    CORE_SRC+=("${DTG_TEXT}.fv_core.res.tile${i}.nc")
    TRACER_SRC+=("${DTG_TEXT}.fv_tracer.res.tile${i}.nc")
    SFC_SRC+=("${DTG_TEXT}.sfc_data.tile${i}.nc")
    ORO_DES+=("${ATMRES}.${OCNRES}_oro_data.tile${i}.nc")
done
CORE_SRC+=("${DTG_TEXT}.fv_core.res.nc")

ORO_SRC=$(printf ',"%s"' "${ORO_SRC[@]}")          && ORO_SRC="${ORO_SRC:2:-1}"
CORE_SRC=$(printf ',"%s"' "${CORE_SRC[@]}")        && CORE_SRC="${CORE_SRC:2:-1}"
TRACER_SRC=$(printf ',"%s"' "${TRACER_SRC[@]}")    && TRACER_SRC="${TRACER_SRC:2:-1}"
SFC_SRC=$(printf ',"%s"' "${SFC_SRC[@]}")          && SFC_SRC="${SFC_SRC:2:-1}"
ORO_DES=$(printf ',"%s"' "${ORO_DES[@]}")          && ORO_DES="${ORO_DES:2:-1}"
export OROG_FILES_INPUT_GRID=${ORO_SRC}
export ATM_CORE_FILES_INPUT=${CORE_SRC}
export ATM_TRACER_FILES_INPUT=${TRACER_SRC}
export SFC_FILES_INPUT=${SFC_SRC}
export OROG_FILES_TARGET_GRID=${ORO_DES}
export COMIN=${dir}
export CDATE=${dtg}
export ocn=${OCNRES:2:3}
export VCOORD_FILE="${HOMEufs}/fix/am/global_hyblev.l128C.txt"
export MOSAIC_FILE_INPUT_GRID="${HOMEufs}/fix/orog/${SRC_ATMRES}/${SRC_ATMRES}_mosaic.nc"
export MOSAIC_FILE_TARGET_GRID="${HOMEufs}/fix/orog/${ATMRES}/${ATMRES}_mosaic.nc"
export OROG_DIR_INPUT_GRID="${HOMEufs}/fix/orog/${SRC_ATMRES}"
export OROG_DIR_TARGET_GRID="${HOMEufs}/fix/orog/${ATMRES}"
export TRACERS_INPUT='"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'
export TRACERS_TARGET=${TRACERS_INPUT}
# from nesmio files?
#export TRACERS_INPUT='"spfh","clwmr","o3mr","icmr","rwmr","snmr","grle"'
#export TRACERS_TARGET='"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'

########################
# modules
module purge
module use ${HOMEufs}/modulefiles
module load build.${m_target}.${compiler}

mkdir -p ${COMIN}/CHGRES
cd ${COMIN}/CHGRES
############
namelist="fort.41"
[[ -f ${namelist} ]] && rm ${namelist}
# namelist options grid
cat <<EOF > ${namelist}

&config
 input_type="restart"
 fix_dir_target_grid="${OROG_DIR_TARGET_GRID}/sfc"
 mosaic_file_target_grid="${MOSAIC_FILE_TARGET_GRID}"
 orog_dir_target_grid="${OROG_DIR_TARGET_GRID}"
 orog_files_target_grid="${OROG_FILES_TARGET_GRID}"
 mosaic_file_input_grid="${MOSAIC_FILE_INPUT_GRID}"
 orog_dir_input_grid="${OROG_DIR_INPUT_GRID}"
 orog_files_input_grid="${OROG_FILES_INPUT_GRID}"
 data_dir_input_grid="${COMIN}"
 atm_core_files_input_grid="${ATM_CORE_FILES_INPUT}"
 atm_tracer_files_input_grid="${ATM_TRACER_FILES_INPUT}"
 vcoord_file_target_grid="${VCOORD_FILE}"
 sfc_files_input_grid="${SFC_FILES_INPUT}"
 cycle_mon=${dtg:4:2}
 cycle_day=${dtg:6:2}
 cycle_hour=${dtg:8:2}
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=.true.
 tracers_input=${TRACERS_INPUT}
 tracers=${TRACERS_TARGET}
/
EOF

${APRUN} -n 12 ${HOMEufs}/exec/chgres_cube
if (( ${?} > 0 )); then
    echo 'chgres_ATM failed'
    exit 1
fi

# move files
mkdir -p ${dir_input_atmos}
for n in $(seq 1 6); do
    mv out.atm.tile${n}.nc ${dir_input_atmos}/gfs_data.tile${n}.nc
    mv out.sfc.tile${n}.nc ${dir_input_atmos}/sfc_data.tile${n}.nc
done
mv gfs_ctrl.nc ${dir_input_atmos}/gfs_ctrl.nc

cd ${dir}
#rm -r CHGRES
#rm ${DTG_TEXT}*.nc
echo "chgres_ATM.sh SUCCESSFUL"
