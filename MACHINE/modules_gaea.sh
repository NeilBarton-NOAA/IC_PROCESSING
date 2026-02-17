#!/bin/sh
set -u
set +x
module purge

module use -a /ncrc/home2/Neil.Barton/TOOLS/modulefiles
module load conda
module use -a /sw/gaea-c6/spack-envs/base/modules/spack/Core
module load nco
module use -a /usw/hpss/modulefiles
module load hsi
set -x
