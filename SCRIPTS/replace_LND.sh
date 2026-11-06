#!/bin/sh
set -u
dtg=${1:-$dtg}
source ${SCRIPT_DIR}/defaults.sh
file_land_states="${dir_land_states}/ufs_land_restart.${dtg:0:4}-${dtg:4:2}-${dtg:6:2}_00-00-00.nc"

PYTHON=$(which python)
${PYTHON} ${SCRIPT_DIR}/LAND_UPDATE/spinup/replace_land_states.py ${dir_input_atmos}/ ${file_land_states} 
[[ $? > 0 ]] && echo "FATAL in replay_land_states.py" && exit 1

echo 'SUCCESSFUL: replace_LAND.sh'
exit 0
