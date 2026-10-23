#!/bin/bash
set -eux

DRUN=${SCRIPT_DIR}/ICE_UPDATE
HEXE=prepare_cice6_restart_driver.py
iconc=0
ithkn=0
hsnow=1
snitd=1
sst=0
fyaml=cice6rest_files_SFSmem.yaml

usage() {
  cat <<EOF
Usage:
  $0 --rdates DATES [options]

Required:
  --rdates YYYYMMDD    restart date(s) to process, may be > than 1 date
  
Optional:
  --iconc 0|1     1 - insert ice conc., 0 - no (default=${iconc})
  --ithkn 0|1     1 - insert ice thickness, 0 - no (default=${ithkn})
  --hsnow 0|1     1 - insert snow depth, 0 - no (default=${hsnow})
  --snitd 0|1     1 - add tracers for ITDrdg, 0 - no (default=${snitd})
  --sst   0|1     1 - update SST in MOM6 restart under sea ice, 0 - no (default=${sst})
  --fyaml FYAML   FYAML file with ice restart input, output, etc. (default=${fyaml})

Example:
  $0 --iconc 1 --ithkn 1 --hsnow 1 --snitd 0 --sst 0 --fyaml ice_config.yaml --rdates 20240101 20240201
EOF
  exit 1
}

# Parse the command-line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --iconc)
      iconc=$2
      shift 2 # Move past the flag and its arg. to the next flag
      ;;
    --ithkn)
      ithkn=$2
      shift 2
      ;;
    --hsnow)
      hsnow=$2
      shift 2
      ;;
    --snitd)
      snitd=$2
      shift 2
      ;;
    --sst)
      sst=$2
      shift 2
      ;;
    --fyaml)
      fyaml=$2
      shift 2
      ;;
    --rdates)
      shift
      RDATES=()
      while [[ $# -gt 0 && $1 != --* ]]; do
        RDATES+=("$1")
        shift
      done
      ;;
    --help)
      usage
      ;;
    *)
    echo "Error: Unrecognized option $1"
    usage
    ;;
  esac
done

if [[ ${#RDATES[@]} -eq 0 ]]; then
  echo "ERROR: no restart dates specified"
  usage
fi

cd $DRUN
pwd

for rdate_in in "${RDATES[@]}"; do
  echo " ==================== "
  echo "Processing ${rdate_in}"
  echo " ==================== "
  python prepare_cice6_restart_driver.py \
      --rdate_in "${rdate_in}" \
      --iconc "${iconc}" \
      --ithkn "${ithkn}" \
      --hsnow "${hsnow}" \
      --snitd "${snitd}" \
      --sst "${sst}" \
      --fyaml "${fyaml}"

done

[[ ${?} > 0 ]] && echo "FATAL ERROR: ICE ICs ARE NOT CORRECTLY REPLACED!" && exit 1

echo "ALL DONE"

# Replace the original ICs with the newly generated ICs
mv ${dir_restart_ice}/tmp/cice_restart.${dtg:0:8}.00.hsnow.snphys.nc ${dir_restart_ice}/${dtg:0:8}.000000.cice_model.res.nc
rm -rf ${dir_restart_ice}/tmp

exit 0
