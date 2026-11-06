#!/bin/sh
set -xu
dtg=${1}
source ${SCRIPT_DIR}/defaults.sh
dir=${dir_restart_ice}
compiler=${chgres_compiler}

########################
HOMEufs=${CODE_DIR}/UFS_UTILS
OCNICEPREP=${HOMEufs}/reg_tests/ocnice_prep/parm
EXEC=${HOMEufs}/exec/oiprep
RT_DIR=${HOMEufs}/reg_tests
source ${RT_DIR}/rt.control 
FIXDIR=${HOMEreg}/cpld_gridgen/baseline_data

########################
if [[ ${OCNRES} == "mx025" ]]; then
    SRCDIMS="360,320"
    DSTDIMS="1440,1080"
elif [[ ${OCNRES} == "mx100" ]]; then
    SRCDIMS="1440,1080"
    DSTDIMS="360,320"
fi

########################
WORKDIR=${dir}/CHGRES
if [[ -d ${WORKDIR} ]]; then
    rm ${WORKDIR}/*
else
    mkdir -p ${WORKDIR} 
fi
cd ${WORKDIR}

ln -sf ${dir}/${DTG_TEXT}.cice_model.res.nc ${WORKDIR}/ice.nc
ln -sf ${OCNICEPREP}/ice.csv ${WORKDIR}

cat << EOF > ocniceprep.nml
&ocniceprep_nml
ftype='ice'
wgtsdir="${FIXDIR}"
griddir="${FIXDIR}"
srcdims=${SRCDIMS}
dstdims=${DSTDIMS}
debug=.true.
/
EOF

########################
# modules
#module purge
[[ "${m_target}" == "ursa" ]] && export sfcio_ver=1.4.2
module use ${HOMEufs}/modulefiles
module load build.${machine}.${compiler}

########################
# run
echo "Running $( basename ${EXEC} ) at ${PWD}"
cp ${EXEC} .
${APRUN} -n 1 ./$( basename ${EXEC} )

if (( ${?} > 0 )); then
    echo 'chgres_ICE failed'
    exit 1
fi
if [[ ! -f ${WORKDIR}/ice.${OCNRES}.nc ]]; then
    echo "FATAL: ${WORKDIR}/ice.${OCNRES}.nc not created"
fi
rm ${dir}/${DTG_TEXT}.cice_model.res.nc
mv ${WORKDIR}/ice.${OCNRES}.nc ${dir}/${DTG_TEXT}.cice_model.res.nc
rm -rf ${WORKDIR}
echo "chgres_ICE.sh SUCCESSFUL"
exit 0
