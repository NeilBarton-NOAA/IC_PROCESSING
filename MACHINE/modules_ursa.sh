#!/bin/sh
set -u
module purge
module load nco
module use -a ~Neil.Barton/TOOLS/modulefiles
module load conda
