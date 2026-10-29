#!/bin/sh
# https://noaa-emcufs-utils.readthedocs.io/en/latest/
# for replay data
set -u
set -x
dtg=${1:-$dtg}
SRC_ATMRES=${2:-"C192"}
SRC_OCNRES=${3:-"mx025"}
MEM=${4:-'000'}
IC_SRC=${5}
NTASKS=${6:-12}
source ${SCRIPT_DIR}/defaults.sh
compiler=${chgres_compiler}
dir=${dir_restart_atmos}
dir=${dir//mem000/mem${MEM}}
dir_input_atmos=${dir_input_atmos//mem000/mem${MEM}}
export DATA=${dir}
cd ${DATA}
echo ${DATA}
########################
# chgres_cube.sh options
export HOMEufs=${CODE_DIR}/UFS_UTILS
set -x
########################
# REPLAY ICs need different ORO files
if [[ ${IC_SRC} == 'REPLAY' ]]; then
    # to fix mismatch in surface files, copy old oro_data files
    #   cp /scratch3/NCEPDEV/global/role.glopara/fix/orog/20230615/C384.mx025/*oro_data* .
    OROG_REPLAY=${HOMEufs}/fix/orog_replay
    if [[ ! -d ${OROG_REPLAY} ]]; then
        OROG_DIR=$(readlink ${HOMEufs}/fix/orog)
        cp -r ${OROG_DIR} ${OROG_REPLAY}
        cp ${OROG_DIR}/../20230615/C384.mx025/*oro_data* ${OROG_REPLAY}/C384
    fi
fi
########################
# needed for IC_SRC == "restart"
ORO_SRC=() && CORE_SRC=() && TRACER_SRC=() && SFC_SRC=() && ORO_DES=()
for i in {1..6}; do
    ORO_SRC+=("${SRC_ATMRES}.${SRC_OCNRES}_oro_data.tile${i}.nc")
    CORE_SRC+=("${DTG_TEXT_MINUS3}.fv_core.res.tile${i}.nc")
    TRACER_SRC+=("${DTG_TEXT_MINUS3}.fv_tracer.res.tile${i}.nc")
    SFC_SRC+=("${DTG_TEXT_MINUS3}.sfc_data.tile${i}.nc")
    ORO_DES+=("${ATMRES}.${OCNRES}_oro_data.tile${i}.nc")
done
CORE_SRC+=("${DTG_TEXT_MINUS3}.fv_core.res.nc")

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
export TRACERS_TARGET='"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'

if [[ ${IC_SRC} == 'REPLAY' ]]; then
    export MOSAIC_FILE_INPUT_GRID="${HOMEufs}/fix/orog_TEST/${SRC_ATMRES}/${SRC_ATMRES}_mosaic.nc"
    export MOSAIC_FILE_TARGET_GRID="${HOMEufs}/fix/orog_TEST/${ATMRES}/${ATMRES}_mosaic.nc"
    export OROG_DIR_INPUT_GRID="${HOMEufs}/fix/orog_TEST/${SRC_ATMRES}"
    export OROG_DIR_TARGET_GRID="${HOMEufs}/fix/orog_TEST/${ATMRES}"
else
    export MOSAIC_FILE_INPUT_GRID="${HOMEufs}/fix/orog/${SRC_ATMRES}/${SRC_ATMRES}_mosaic.nc"
    export MOSAIC_FILE_TARGET_GRID="${HOMEufs}/fix/orog/${ATMRES}/${ATMRES}_mosaic.nc"
    export OROG_DIR_INPUT_GRID="${HOMEufs}/fix/orog/${SRC_ATMRES}"
    export OROG_DIR_TARGET_GRID="${HOMEufs}/fix/orog/${ATMRES}"
fi
if [[ ${IC_SRC} == *"CPC"* ]]; then
    export TRACERS_INPUT='"spfh","clwmr","o3mr","icmr","rwmr","snmr","grle"'
    export INPUT_TYPE="gaussian_nemsio"
    export ATM_FILES_INPUT="sfg_${dtg}_fhr06_mem012"
    export SFC_FILES_INPUT="bfg_${dtg}_fhr06_mem012"
    export CONVERT_NST=".false."
    EXEC=chgres_cube_nesmio
else
    export TRACERS_INPUT='"sphum","liq_wat","o3mr","ice_wat","rainwat","snowwat","graupel"'
    export INPUT_TYPE="restart"
    export ATM_FILES_INPUT="None"
    export CONVERT_NST=".true."
    EXEC=chgres_cube_restart
fi

########################
# modules
[[ "${m_target}" == "ursa" ]] && export sfcio_ver=1.4.2
if [[ ${m_target} == "wcoss2" ]]; then
    module reset
else
    module purge
fi
module use ${HOMEufs}/modulefiles
module load build.${m_target}.${compiler}
module list

mkdir -p ${COMIN}/CHGRES
cd ${COMIN}/CHGRES
############
namelist="fort.41"
[[ -f ${namelist} ]] && rm ${namelist}
# namelist options grid
cat <<EOF > ${namelist}
&config
 input_type="${INPUT_TYPE}"
 fix_dir_target_grid="${OROG_DIR_TARGET_GRID}/sfc"
 mosaic_file_target_grid="${MOSAIC_FILE_TARGET_GRID}"
 orog_dir_target_grid="${OROG_DIR_TARGET_GRID}"
 orog_files_target_grid="${OROG_FILES_TARGET_GRID}"
 vcoord_file_target_grid="${VCOORD_FILE}"
 data_dir_input_grid="${COMIN}"
 sfc_files_input_grid="${SFC_FILES_INPUT}"
 cycle_mon=${dtg:4:2}
 cycle_day=${dtg:6:2}
 cycle_hour=${dtg:8:2}
 convert_atm=.true.
 convert_sfc=.true.
 convert_nst=${CONVERT_NST}
 tracers_input=${TRACERS_INPUT}
 tracers=${TRACERS_TARGET}
EOF
if [[ ${INPUT_TYPE} == 'restart' ]]; then
cat <<EOF >> ${namelist}
 mosaic_file_input_grid="${MOSAIC_FILE_INPUT_GRID}"
 orog_dir_input_grid="${OROG_DIR_INPUT_GRID}"
 orog_files_input_grid="${OROG_FILES_INPUT_GRID}"
 atm_core_files_input_grid="${ATM_CORE_FILES_INPUT}"
 atm_tracer_files_input_grid="${ATM_TRACER_FILES_INPUT}"
/
EOF
else
cat <<EOF >> ${namelist}
 atm_files_input_grid="${ATM_FILES_INPUT}"
/
EOF
fi
${APRUN} -n ${NTASKS} ${HOMEufs}/exec/${EXEC}
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
rm -r CHGRES
if [[ ${INPUT_TYPE} == 'restart' ]]; then
    rm ${DTG_TEXT_MINUS3}*.nc
else
    rm ${ATM_FILES_INPUT} ${SFC_FILES_INPUT}
fi
cd ${dir_restart_atmos}/../../
rm -rf atmos

echo "chgres_ATM.sh SUCCESSFUL"
