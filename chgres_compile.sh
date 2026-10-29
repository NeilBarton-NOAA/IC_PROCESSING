#!/bin/sh
set -u
IC_SRC=${1:-'GFS'}
# compiles chgres program
SCRIPT_DIR=$(dirname "$0")/SCRIPTS
dtg=2017100100 # dummy variable
source ${SCRIPT_DIR}/defaults.sh
CODE=https://github.com/ufs-community/UFS_UTILS.git
HASH=develop
export compiler=${chgres_compiler}
export target=${m_target}

mkdir -p ${CODE_DIR}
cd ${CODE_DIR}
git clone ${CODE}
cd ${CODE_DIR}/UFS_UTILS
git checkout ${HASH}

# fix files
cd ${CODE_DIR}/UFS_UTILS/fix
bash link_fixdirs.sh emc ${m_target}

# edit build
if [[ ${m_target} == "gaea" ]]; then
    build_file=${CODE_DIR}/UFS_UTILS/modulefiles/build.${target}c6.${compiler}.lua
else
    build_file=${CODE_DIR}/UFS_UTILS/modulefiles/build.${target}.${compiler}.lua
fi
sed -i 's/--sfcio_ver/sfcio_ver/' ${build_file}
sed -i 's/--load(pathJoin("sfcio"/load(pathJoin("sfcio"/' ${build_file}
mkdir -p ${CODE_DIR} && cd ${CODE_DIR}

# compile
cd ${CODE_DIR}/UFS_UTILS

## build for GFS/RESART ICs
[[ "${m_target}" == "ursa" ]] && export sfcio_ver=1.4.2
bash build_all.sh
if (( ${?} != 0 )); then
    echo 'COMPILE failed'
    exit 1
fi
cp ${CODE_DIR}/UFS_UTILS/exec/chgres_cube ${CODE_DIR}/UFS_UTILS/exec/chgres_cube_restart

# build for CPC ICs nesmio
export CMAKE_OPTS="-DCHGRES_ALL=ON"
bash build_all.sh
if (( ${?} != 0 )); then
    echo 'COMPILE failed'
    exit 1
fi
cp ${CODE_DIR}/UFS_UTILS/exec/chgres_cube ${CODE_DIR}/UFS_UTILS/exec/chgres_cube_nesmio

