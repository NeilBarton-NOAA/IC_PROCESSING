"""
  After correcting CICE6 initial fields (ice concentration)
  Updated MOM6 surface T:
  under sea ice T = Tfreeze(S)
  for aice < 1: weighted update 

   Need MOM restart with surf T and CICE6 restart with corrected aice
   to use for updating SST under the sea ice

   author: Dmitry Dukhovskoy NOAA NWS OMD, May 2026
"""
import os
import numpy as np
import sys
import importlib
import matplotlib
import xarray
import matplotlib.colors as colors
import yaml
import argparse
from pathlib import Path
import sys
                      
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from MyPython.mod_swstate import freezing_temp_linear, sw_dens0
from MyPython.mod_mom6 import read_mom6grid

CONFIG_DIR = ROOT / "config"
default_yaml = next(CONFIG_DIR.glob("*.yaml")).name

print("\n ======= MOM6 correcting SST  ==========")

yrR = mmR = ddR = hrR = None
yrN = mmN = ddN = hrN = None

# Current code assumes 1 MOM restart file
# No date / time change 
# TODO: add logic for N MOM restart files: MOM_*.res.nc, MOM_*.res_1.nc, ... 
# YAML is used if input file name or vartmp is None
#
# flice - CICE6 restart file with updated ice conc. fields
#              where SST will be updated
#         provided by the driver script or (if not) read from YAML

vartmp = 'Temp'      # variable name in the restart file

parser = argparse.ArgumentParser()
parser.add_argument("--flmom_in", help=f"MOM6 rest file name, original", type=str)
parser.add_argument("--flice", help="CICE6 restart file to use for MOM6 surf T", type=str)
parser.add_argument("--fyaml",
    help=f"YAML with local directories, filenames, restart dates, default={default_yaml}, inside ./config",
    default=default_yaml,
    type=str)
args = parser.parse_args()

flmom_in   = args.flmom_in if args.flmom_in else None 
flrst_cice = None if args.flice in (None, "None", "") else args.flice
fyaml      = CONFIG_DIR / args.fyaml

print(f"Reading YAML with restart info: {fyaml}\n")
with open(fyaml) as ff:
  config_rest = yaml.safe_load(ff)

def expand_env_vars(d):
    for k, v in d.items():
        if isinstance(v, dict):
            expand_env_vars(v)
        elif isinstance(v, str):
            d[k] = os.path.expandvars(v)

expand_env_vars(config_rest)

if vartmp is None:
  vartmp = config_rest["mom6_names"]["varnm_temp"]
  assert vartmp is not None, "vartmp is not provided"

if flmom_in is None:
  flmom_in = config_rest["mom6_names"]["flmom_in"]

assert flmom_in is not None, f"flmom_in has to be provided either --flmom_in or in {fyaml}"

# MOM restart path: input and output restart dir
pthrest_in  = config_rest["mom6_paths"]["pth_in"]
pthrest_out = config_rest["mom6_paths"]["pth_out"]
if pthrest_out is None:
  pthrest_out = pthrest_in 

assert pthrest_in is not None, f"pthrest_in is not specified in {fyaml}"
assert pthrest_out is not None, f"pthrest_out is None"

os.makedirs(pthrst_out, exist_ok=True)

# Get MOM6 grid
pthgrid = config_rest["mom6_paths"]["pth_grid"]
dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

with xarray.open_dataset(dftopo_mom) as dstopo:
  HH = dstopo['depth'].data.squeeze()

HH = np.where(HH < 1.e-20, np.nan, HH)
HH = -HH
HH = np.where(np.isnan(HH), 1., HH)
jdm, idm = HH.shape 

LON, LAT = read_mom6grid(dfgrid_mom, grdpnt='hgrid')

# CICE parameters:
puny      = 1.e-11

# CICE6 restart file:
pthcice_in = config_rest["mom6_paths"]["pth_ice"]
# If not provided, use specified in YAMLS
if flrst_cice is None:
  print("Reading flrst_cice")
  flrst_cice = config_rest["mom6_names"]["flrst_ice"]

assert flrst_cice is not None, f"flrst_cice is None, check {fyaml} or use --flice"
assert pthcice_in is not None, f"pthcice_in is None, check {fyaml} mom6_names"
 
dflcice = os.path.join(pthcice_in, flrst_cice)
print(f"Reading CICE6 restart: {dflcice}")
with xarray.open_dataset(dflcice) as ds_ice:
  aicen = ds_ice["aicen"].data  # partial ice area by cats

# Aggregated ice partial area:
aice = np.sum(aicen, axis=0).squeeze()

# MOM input restart:
dflmom_in = os.path.join(pthrest_in, flmom_in) 
print(f"Reading restart: {dflmom_in}")
ds_in = xarray.open_dataset(dflmom_in, decode_times=False)
ds_out = ds_in.copy(deep=True)

# 3D fields:
varlr  = 'eta'   # interface depth
varsal = 'Salt'  # salinity

temp_in  = ds_in[vartmp].isel(Time=0).data  # 3D pot temp
lrthk_in = ds_in[varlr].isel(Time=0).data  # interf. depths, lr0 = surf interf
salt_in   = ds_in[varsal].isel(Time=0).data  

ds_in.close()

# Check signs:
assert np.nanmin(HH) < 0, f"Bottom topo has to be negative"
assert np.nanmin(lrthk_in[1,:,:]) < 0, f"Layer interf depths expected <0"

# Compute freezing T:
S2d = salt_in[0,:,:].squeeze()
S2d[HH >= 0] = np.nan
Tfrz = freezing_temp_linear(S2d,0)

# Update surf SST in MOM6:
T2d = temp_in[0,:,:].squeeze()

T2d_new = T2d*(1. - aice) + aice * Tfrz
T2d_new = np.where(aice < puny, T2d, T2d_new)

T3d_new = temp_in.copy()
T3d_new[0,:,:] = T2d_new
S3d_new = salt_in.copy()

def drho_lrs(ilr, S3d, T3d, aice, puny, HH, lrthk_in):
  rho1 = sw_dens0(S3d[ilr,:,:], T3d[ilr,:,:])
  rho2 = sw_dens0(S3d[ilr+1,:,:], T3d[ilr+1,:,:])
  drho = rho1 - rho2   
  # Mak not ice regions:
  drho[aice < puny] = np.nan
  # Mask regions where layer depth is shallower than local bottom
  # HH < 0 for ocean grid points
  zzLr = lrthk_in[ilr,:,:]  # upper interface
  drho[zzLr < HH] = np.nan 

  return drho

def mix_layers(ilr, lrthk_in, HH, T3d_new, S3d_new, JJ, II):
  # Interface depths:
  zlr1 = lrthk_in[ilr,:,:]
  zlr2 = lrthk_in[ilr+1,:,:]
  zlr3 = lrthk_in[ilr+2,:,:]
  
  # mask out below-bottom vanished layers
  # last interf = bottom depth
  zlr1 = np.where(zlr1 < HH, HH, zlr1)
  zlr1[HH >= 0] = np.nan
  zlr2 = np.where(zlr2 < HH, HH, zlr2)
  zlr2[HH >= 0] = np.nan
  zlr3 = np.where(zlr3 < HH, HH, zlr3)
  zlr3[HH >= 0] = np.nan

  # Lyaer thickness:
  dz1 = np.abs(zlr2 - zlr1)
  dz2 = np.abs(zlr3 - zlr2)
 
  # mix two layers where drho > eps:
  T1 = T3d_new[ilr,:,:]
  T2 = T3d_new[ilr+1,:,:]
  S1 = S3d_new[ilr,:,:]
  S2 = S3d_new[ilr+1,:,:]
  den = dz1 + dz2
  den = np.where(den == 0, np.nan, den)
  tmix =  (dz1*T1 + dz2*T2) / den
  smix =  (dz1*S1 + dz2*S2) / den

  T3d_new[ilr,JJ,II]   = tmix[JJ,II]
  T3d_new[ilr+1,JJ,II] = tmix[JJ,II]
  S3d_new[ilr,JJ,II]   = smix[JJ,II]
  S3d_new[ilr+1,JJ,II] = smix[JJ,II]

  return T3d_new, S3d_new
    

# Next: check water column stability in the upper ocean:
# Potential density ref = 0, assumed that only upper ocean needs adjustment
# Note that in the original fields, some layers can be unstable which seems to be ok for MOM6
# Therefore, do not mix all whatercolumn but only the upper few layers if needed
eps_drho = 1e-2    # threshold for delta rho
mixlrs = 3     # the upper layers to mix if unstable after Tfrz 
nlrs = salt_in.shape[0]
ilr = 0
max_iter = 50
print(f"Adjusting unstable layers under sea ice, max_iter={max_iter}")
for iter in range(max_iter):
  unstb = False

  print(f"iter = {iter+1}")
  for ilr in range(mixlrs):
    drho = drho_lrs(ilr, S3d_new, T3d_new, aice, puny, HH, lrthk_in)
    print(f"  Lr={ilr+1:02d}, max drho={np.nanmax(drho):.6f}")
 
    JJ, II = np.where(drho > eps_drho)
    if JJ.size:
      unstb = True

      # mix layers ilr and ilr+1
      T3d_new, S3d_new = mix_layers(ilr, lrthk_in, HH, T3d_new, S3d_new, JJ, II)      

  if not unstb:
    print(f"Converged in {iter+1} iterations")
    break

if unstb:
  raise RuntimeError("Instability adjustment did not converge after {iter+1} iterations")

# Check for NaNs:
if np.isnan(T3d_new).any():
  print("T3d Converting NaNs --> 0")
  T3d_new = np.nan_to_num(T3d_new, nan=0.0)

if np.isnan(S3d_new).any():
  print("S3d Converting NaNs --> 0")
  S3d_new = np.nan_to_num(S3d_new, nan=0.0)


T3d_new = np.expand_dims(T3d_new, axis=0)
S3d_new = np.expand_dims(S3d_new, axis=0)
ds_out[vartmp].values[:] = T3d_new
ds_out[varsal].values[:] = S3d_new


ds_out.attrs.update({
    "title": f"MOM6 restart with adjusted SST under ice to Tfrz weighted by aice {flrst_cice}",
    "source": "correct_surfT_mom6_restart.py",
})

# Save:flmom_in
# Construct output file name if not provided
# use cice_restart used for aice
# assumed *.res.[sfx].nc  where sfx is iconc, hsnow, etc
# if missing - control run
flrst_out = config_rest["mom6_names"]["flmom_out"]
if flrst_out is None:
  base = os.path.splitext(flmom_in)[0]
  cice_parts = flrst_cice.split(".")
  cice_sfx = cice_parts[-2]
  if cice_sfx == 'res':
    cice_sfx = 'cntrl'

  flrst_out = f"{base}.{cice_sfx}.Tfrz.nc"

dflrst_out = os.path.join(pthrest_out, flrst_out)
print(f"Saving MOM6 restart --> {dflrst_out}")
ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
ds_out.close()


