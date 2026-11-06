#!/bin/sh
set -xu
dtg=${1}
source ${SCRIPT_DIR}/defaults.sh
dir=${dir_restart_ocean}
compiler=${chgres_compiler}

########################
HOMEufs=${CODE_DIR}/UFS_UTILS
OCNICEPREP=${HOMEufs}/reg_tests/ocnice_prep/parm
EXEC=${HOMEufs}/exec/oiprep
RT_DIR=${HOMEufs}/reg_tests
source ${RT_DIR}/rt.control 
FIXDIR=${HOMEreg}/cpld_gridgen/baseline_data

########################
WORKDIR=${dir}/CHGRES
if [[ -d ${WORKDIR} ]]; then
    rm ${WORKDIR}/*
else
    mkdir -p ${WORKDIR} 
fi
cd ${WORKDIR}

if [[ ${OCNRES} == "mx025" ]]; then
    SRCDIMS="360,320"
    DSTDIMS="1440,1080"
elif [[ ${OCNRES} == "mx100" ]]; then
    SRCDIMS="1440,1080"
    DSTDIMS="360,320"
fi

if [[ ! -f ${WORKDIR}/ocean.nc ]]; then
    echo 'Creating ${WORKDIR}/ocean.nc file'
    if [[ ${IC_SRC} == *"CPC"* ]]; then
        ncks -v Temp,Salt,h,sfc,u,v ${dir}/${DTG_TEXT}.MOM.res.nc ${WORKDIR}/ocean.nc
    elif [[ ${IC_SRC} == REPLAY ]]; then
        ncks -v Temp,Salt,h,u ${dir}/${DTG_TEXT}.MOM.res.nc ${WORKDIR}/ocean.nc
        ncks -v v,sfc -A ${dir}/${DTG_TEXT}.MOM.res_1.nc ${WORKDIR}/ocean.nc
    fi
fi
ln -sf ${OCNICEPREP}/ocean.csv ${WORKDIR}

cat << EOF > ocniceprep.nml
&ocniceprep_nml
ftype='ocean'
wgtsdir="${FIXDIR}"
griddir="${FIXDIR}"
srcdims=${SRCDIMS}
dstdims=${DSTDIMS}
debug=.true.
/
EOF

########################
# modules
[[ "${m_target}" == "ursa" ]] && export sfcio_ver=1.4.2
module purge
module use ${HOMEufs}/modulefiles
module load build.${machine}.${compiler}

echo "Running $( basename ${EXEC} ) at ${PWD}"
cp ${EXEC} .
${APRUN} -n 1 ./$( basename ${EXEC} )
if (( ${?} > 0 )); then
    echo 'chgres_OCN failed'
    exit 1
fi
if [[ ! -f ${WORKDIR}/ocean.${OCNRES}.nc ]]; then
    echo "FATAL: ${WORKDIR}/ocean.${OCNRES}.nc not created"
    exit 1
fi

mv ${WORKDIR}/ocean.${OCNRES}.nc ${dir}/${DTG_TEXT}.MOM.res.nc
rm -rf ${WORKDIR}
if [[ ${OCNRES} == "mx100" ]]; then
    rm ${dir}/${DTG_TEXT}.MOM.res_*.nc
fi
echo "chgres_OCN.sh SUCCESSFUL"
exit 0
