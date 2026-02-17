#!/bin/sh
machine=$(uname -n)
echo "${SCRIPT_DIR}/MACHINE/modules_${machine}.sh"
source ${SCRIPT_DIR}/MACHINE/modules_${machine}.sh
