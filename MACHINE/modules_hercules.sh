#!/bin/sh
set -u
set +x
module purge
module use -a $HOME/TOOLS/modulefiles
module load conda

module load nco
set -x
