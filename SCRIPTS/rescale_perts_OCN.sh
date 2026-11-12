#!/bin/bash
# ==============================================================================
# Driver Script, Re-center and Rescale CPC Ocean Perturbations
# ==============================================================================
set -euo pipefail
# Command-line input for IC date (optional, defaults to 20260201 if empty)
INPUT_DATE="${1:-2026020100}"
# Config Options
export YYYYMMDD="${INPUT_DATE:0:8}"     # Target date (defaults to 20260201 if $1 empty)
export HH="00"                          # Cycle hour (00, 06, 12, 18)
##################
# Ensemble Options 
# Set member bounds to process (e.g., 1 to 10, 1 to 30, or 11 to 20)
export MEM_START=1
export MEM_END=10
export MEM_FORMAT="%03d"                # Zero-padded 3-digit format (mem001, mem002...)
# Calculate total number of members dynamically
export NUM_MEMBERS=$(( MEM_END - MEM_START + 1 ))

##################
# Path & Directory Definitions
# Raw original CPC perturbations root directory
export OCEAN_PERT_ORIGIN=${OCEAN_PERT_ORIGIN:-"/scratch4/NCEPDEV/global/Yangxing.Zheng/ICs/CPC_NRT/C192mx025"}
# Target cycle directory path
export CYCLE_ORIGIN_DIR="${OCEAN_PERT_ORIGIN}/sfs.${YYYYMMDD}/${HH}"
# Operational Working Directory (Scratch Space: work/YYYYMMDD/HH)
export WORK_DIR="${CYCLE_ORIGIN_DIR}/WORK"
export OCEAN_PERT_RESCALED=${CYCLE_ORIGIN_DIR}/RESCALED

# Dynamically extract the exact perturbation filename from the first member directory
FIRST_MEM_DIR=$(printf "${CYCLE_ORIGIN_DIR}/mem${MEM_FORMAT}/analysis/ocean" "${MEM_START}")
export PERT_FILENAME=$(basename $(ls ${FIRST_MEM_DIR}/*.mom6_perturbation.nc 2>/dev/null | head -n 1))
# Path for saving computed ensemble mean file inside working directory
export ENSEMBLE_MEAN_FILE="${WORK_DIR}/ensemble_mean.nc"
# Fixed Scaling Directory & Profile CSV
FIXED_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)/ocn_perts_scale_factors"
export SCALE_CSV="${FIXED_DIR}/scale.csv"

##################
# Constraints
export MAX_SCALE_RATIO=3.0              # Maximum scaling factor cap (prevents numerical instability)
export NUM_WORKERS="${NUM_MEMBERS}"     # One worker process per ensemble member
export NETCDF_COMPRESS="YES"            # Set "YES" for zlib compression, "NO" for max speed
export KEEP_ORIG=${KEEP_ORIG:-"NO"}    # KEEP original rescaled perturbations
export CLEAN_WORK_DIR="YES"              # Set "YES" to purge temporary work space on success

##################
# Print Workflow Environment Check
echo "======================================================================"
echo "      Driver Script, Re-center and Rescale CPC Ocean Perturbations    "
echo "======================================================================"
echo "  Target IC Date       : ${YYYYMMDD}"
echo "  Cycle Hour           : ${HH}Z"
echo "  Ensemble Range       : mem$(printf "${MEM_FORMAT}" "${MEM_START}") to mem$(printf "${MEM_FORMAT}" "${MEM_END}") (${NUM_MEMBERS} members)"
echo "  Perturbation File    : ${PERT_FILENAME:-[Not Found]}"
echo "  Source Directory     : ${CYCLE_ORIGIN_DIR}"
echo "  Scaling Profile CSV  : ${SCALE_CSV}"
echo "  Max Scaling Cap      : ${MAX_SCALE_RATIO}x"
echo "  Work Space           : ${WORK_DIR}"
echo "  Ensemble Mean Output : ${ENSEMBLE_MEAN_FILE}"
echo "  Final Output Dir     : ${OCEAN_PERT_RESCALED}/sfs.${YYYYMMDD}/${HH}"
echo "======================================================================"

##################
# 3. Strict Pre-flight Checks
if [[ ! -d "${CYCLE_ORIGIN_DIR}" ]]; then
    echo "[FATAL ERROR] Target cycle origin directory does not exist:"
    echo "              ${CYCLE_ORIGIN_DIR}"
    exit 1
fi
if [[ -z "${PERT_FILENAME:-}" ]]; then
    echo "[FATAL ERROR] Could not locate any .mom6_perturbation.nc file in:"
    echo "              ${CYCLE_ORIGIN_DIR}"
    exit 1
else
    echo "✓ Found source perturbation file: ${PERT_FILENAME}"
fi
if [[ ! -f "${SCALE_CSV}" ]]; then
    echo "[FATAL ERROR] Required fixed scale CSV file is missing:"
    echo "              ${SCALE_CSV}"
    exit 1
else
    echo "✓ Found fixed scale profile: ${SCALE_CSV}"
fi

##################
# Create Working and Output Directories
mkdir -p "${WORK_DIR}"
mkdir -p "${OCEAN_PERT_RESCALED}"
echo "✓ Directories initialized successfully."

##################
# Compute Ensemble Mean for t_pert and s_pert
echo ""
PYTHON_MEAN_SCRIPT="${SCRIPT_DIR}/pytools/compute_ensemble_mean.py"
if [[ ! -f "${PYTHON_MEAN_SCRIPT}" ]]; then
    echo "[FATAL ERROR] Python mean calculation script missing: ${PYTHON_MEAN_SCRIPT}"
    exit 1
fi
python "${PYTHON_MEAN_SCRIPT}" \
    "${CYCLE_ORIGIN_DIR}" \
    "${MEM_START}" \
    "${MEM_END}" \
    "${PERT_FILENAME}" \
    "${ENSEMBLE_MEAN_FILE}"

##################
# Parallel Re-centering (Subtract Ensemble Mean)
echo ""
PYTHON_RECENTER_SCRIPT="${SCRIPT_DIR}/pytools/recenter_parallel.py"
if [[ ! -f "${PYTHON_RECENTER_SCRIPT}" ]]; then
    echo "[FATAL ERROR] Python recentering script missing: ${PYTHON_RECENTER_SCRIPT}"
    exit 1
fi
python "${PYTHON_RECENTER_SCRIPT}" \
    "${CYCLE_ORIGIN_DIR}" \
    "${WORK_DIR}" \
    "${MEM_START}" \
    "${MEM_END}" \
    "${PERT_FILENAME}" \
    "${ENSEMBLE_MEAN_FILE}" \
    "${NUM_WORKERS}"

##################
# Parallel Rescaling (Apply Profile Scaling with Depth Interpolation & Capping)
echo ""
PYTHON_RESCALE_SCRIPT="${SCRIPT_DIR}/pytools/rescale_parallel.py"
if [[ ! -f "${PYTHON_RESCALE_SCRIPT}" ]]; then
    echo "[FATAL ERROR] Python rescaling script missing: ${PYTHON_RESCALE_SCRIPT}"
    exit 1
fi
python "${PYTHON_RESCALE_SCRIPT}" \
    "${WORK_DIR}" \
    "${OCEAN_PERT_RESCALED}" \
    "${YYYYMMDD}" \
    "${HH}" \
    "${MEM_START}" \
    "${MEM_END}" \
    "${PERT_FILENAME}" \
    "${SCALE_CSV}" \
    "${MAX_SCALE_RATIO}" \
    "${NUM_WORKERS}"

##################
# Move and re-write files back to original directory
echo ""
echo "======================================================================"
echo " COPY FILES TO RUN DIRECTORY"
echo "======================================================================"
DES_DIR="${OCEAN_PERT_ORIGIN}/sfs.${YYYYMMDD}/${HH}"
SRC_DIR="${OCEAN_PERT_RESCALED}"
KEEP_ORIG=NO
if [[ ${KEEP_ORIG} == 'YES' ]]; then
    files=$(find ${DES_DIR} -name *mom6_perturbation.nc)
    for f in ${files}; do
        cp ${f} ${f}_ORIG
    done
fi
files=$(find ${SRC_DIR} -name "*mom6_perturbation.nc")
for f in ${files}; do
    printf 'cp \t%s\n\t%s/mem%s\n' "$f" "$DES_DIR" "${f#*mem}"
    cp "${f}" "${DES_DIR}/mem${f#*mem}"
done

##################
# CLEANUP (Optional)
if [[ "${CLEAN_WORK_DIR}" == "YES" ]]; then
    echo ""
    echo "Purging temporary work directories: "
    rm -rf "${WORK_DIR}"
    rm -rf "${OCEAN_PERT_RESCALED}"
fi

echo ""
echo "======================================================================"
echo "✓ ALL STEPS COMPLETED SUCCESSFULLY FOR DATE: ${YYYYMMDD}"
echo "======================================================================"
