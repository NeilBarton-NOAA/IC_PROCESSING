#!/usr/bin/env python3 
########################
#  Neil P. Barton (NOAA-EMC), 2022-10-27
########################
import argparse
import numpy as np
import shutil
import xarray as xr

################################################
# parse arguments
parser = argparse.ArgumentParser( description = "edit sfce file")
parser.add_argument('-f', '--file', action = 'store', nargs = 1, \
        help="sfc_data.tile?.nc file that needs to be edited")
args = parser.parse_args()
f = args.file[0]
print('Surface File:', f)
dat = xr.open_dataset(f)

################################################
# change canopy value
var = 'canopy' 
data = dat[var].values
data[data > 0.5] = 0.5
dat[var] = (dat[var].dims, data)
dat.assign_attrs({'File Edited' : 'canopy values reset to 0.5'})

################################################
# remove NaNs
VARS = set(dat.variables.keys())
for var in VARS:
    data = dat[var].values
    data[np.isnan(data) == True] = 1e30
    dat[var] = (dat[var].dims, data)
    
################################################
# copy file and save new file
shutil.move(f, f + '_old')
dat.to_netcdf(f, mode = 'w')
print('SAVED: ', f)

