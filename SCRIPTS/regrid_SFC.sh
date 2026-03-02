#!/bin/sh
# https://noaa-emcufs-utils.readthedocs.io/en/latest/
set -u
dtg=${1}
SRC_ATMRES=${2:-'C1152'}
MEM=${3:-'000'}
DES_OCNRES=${ONCRES:-'mx025'}
DES_ATMRES=${ATMRES:-'C192'}
source ${SCRIPT_DIR}/defaults.sh
compiler=${chgres_compiler}
dir=${dir_inc_atmos}
dir=${dir//mem000/mem${MEM}}
export HOMEufs=${CODE_DIR}/UFS_UTILS
FIXorog=${HOMEufs}/fix/orog
LONB_CASE_IN=$((4*${SRC_ATMRES:1}))
LATB_CASE_IN=$((2*${SRC_ATMRES:1}))

####################################
DATA=${dir}/REGRID 
[[ -d ${DATA} ]] && rm -r ${DATA}
mkdir -p ${DATA} && cd ${DATA}
echo "DATA dir ${DATA}"
set -x
if [[ ${MEM} == "000" ]]; then
    NAME="ensmean_increment"
    new_dir="${dir/mem000\//mem000\/ensstat\/}" 
    cp ${new_dir}/*increment.sfc.i*.nc .
else
    NAME="increment"
    cp ${dir}/*increment.sfc.i*.nc .
fi
####################################
# set up fix files
ln -sf "${FIXorog}/${SRC_ATMRES}/gaussian.${LONB_CASE_IN}.${LATB_CASE_IN}.nc" "${DATA}/gaussian_scrip.nc"
ln -sf "${FIXorog}/${DES_ATMRES}/${DES_ATMRES}_mosaic.nc" "${DATA}/${DES_ATMRES}_mosaic.nc"
ntiles=6
for n in $(seq 1 $ntiles); do
    ln -sf ${FIXorog}/${DES_ATMRES}/sfc/${DES_ATMRES}.${DES_OCNRES}.vegetation_type.tile${n}.nc  ${DATA}/vegetation_type.tile${n}.nc
    ln -sf ${FIXorog}/${DES_ATMRES}/${DES_ATMRES}_grid.tile${n}.nc ${DATA}/${DES_ATMRES}_grid.tile${n}.nc
done

####################################
# namelist
cat << EOF > regrid.nml
 &config
  n_vars=4,
  variable_list="soilt1_inc", "soilt2_inc", "slc1_inc", "slc2_inc",
  missing_value=0.,
  time_list=3,6,9,
  add_time_dim=.true.,
  extrap_levs=2,
  nmem_ens=1,
 /
 &input
  gridtype="gau_inc",
  ires=${LONB_CASE_IN},
  jres=${LATB_CASE_IN},
  fname="sfs.t00z.${NAME}.sfc.i",
  dir="./",
  fname_coord="gaussian_scrip.nc",
  dir_coord="./"
 /

 &output
  gridtype="fv3_rst",
  ires=${DES_ATMRES:1},
  jres=${DES_ATMRES:1},
  fname="increment.sfc",
  dir="./",
  fname_mask="vegetation_type"
  dir_mask="./"
  dir_coord="./",
 /
EOF

# modules
module purge
module use ${HOMEufs}/modulefiles
module load build.${m_target}.${compiler}

${APRUN} -n 6 ${HOMEufs}/exec/regridStates.x
if (( ${?} > 0 )); then
    echo 'regrid_SFC failed'
    exit 1
fi

files=$(ls increment.sfc.mem*nc)
for f in ${files}; do
    echo "mv ${DATA}/${f} ${dir}/${f//.mem1/}"
    mv ${DATA}/${f} ${dir}/${f//.mem1/}
done
echo "regrid_SFC.sh SUCCESSFUL"
