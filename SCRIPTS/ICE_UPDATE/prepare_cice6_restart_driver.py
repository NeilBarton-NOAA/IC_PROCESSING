"""
  Driver to prepare CICE6 restart for specified
  setup of initial conditions

  Can use any restart adding missing fields to it, e.g.:
  edit YAML

  Specify restart input date/ hour, options:
  (1) use command-line argument to override YAML or if YAML has ~ for rdate_in and / or rhr_in
  (2) provide rdate_in, rhr_in in YAML (and do not use command-line arg.)
  (3) do not use command-line arg., have ~ in YAML, have to provide correct flrst_in in YAML
       then date / time will be deterimed from netcdf attributes (if exist)



  rdate_out:  requested restart date in the new restart
  rdate_in, rhr_in: can be blank, the date/time will be derived from the input restart name assuming one of the formats:
                    YYYYMMDD.XX[XXXX].restart_name.*.nc
                    some_sfx.YYYYMMDD.XX[XXXX].*.nc

  pth_in, pth_out: specify input/output restart dirctories

  flrst_in: intput restart file name
  flrst_out:  leave blank, output restart name will be created based on the flrst_tmp (template) and restart time
               by adding fields being corrected, e.g. cice_restart.20240701.00.iconc_ithkn.nc

  Example:
  create restart with corrected hsnow:
  in YAML file:
  rdate_in:  ~          # input restart date
  rhr_in:    ~          # input restart hour
  rdate_out: 20240701   # output restart date

  flrst_in: "20240701.064800.cice_restart.nc"
  flrst_out: ~
  flrst_tmp: "cice_restart.YYYYMMDD.HH"

  run prepare_cice6_restart_driver.py --iconc 0 --ithkn 0 --hsnow 1 --snitd 0 --regn global

  To create restart with iconc + ithkn + hsnow
  2 options:
  (1) use ice restart from flrst_in (with already inserted iconc and ithkn) and will add
    hsnow on top of this
  flrst_in: "cice_restart.20240701.00.iconc_ithkn.nc"   <--- use restart with already added iconc and ithkn
  flrst_out: ~
  flrst_tmp: "cice_restart.YYYYMMDD.HH"

  run prepare_cice6_restart_driver.py --iconc 0 --ithkn 0 --hsnow 1 --snitd 0 --regn global

  (2) create from original restart, note that intermediate restart files will also be created and saved
  flrst_in: "20240701.064800.cice_restart.nc"   <--- use original restart
  flrst_out: ~
  flrst_tmp: "cice_restart.YYYYMMDD.HH"

  run prepare_cice6_restart_driver.py --iconc 1 --ithkn 1 --hsnow 1 --snitd 0 --regn global


  Note that for snitd = 1, use ice_in with added options for snow redistribution, snow physics

  author: Dmitry Dukhovskoy NOAA NWS OMD, May 2026
"""
import os
import numpy as np
import importlib
import sys
import xarray
import yaml
import argparse
import subprocess
from pathlib import Path

#import MyPython.mod_cice6_utils as mc6
#importlib.reload(mc6)

from MyPython.mod_cice6_utils import change_base_template, flname_replace_date
from MyPython.mod_time import addhours_date

ROOT = Path(__file__).resolve().parent
CONFIG_DIR = ROOT / "config"
CODE_DIR = ROOT / "insertion_cice6"

default_yaml = next(CONFIG_DIR.glob("*.yaml")).name


rest_hr   = 0
rhr_out   = rest_hr
flrst_in  = None
flrst_out = None
#fyaml = 'cice6rest_files_SFSmem.yaml'  # to run Neils restart for SFS ensmb GFS / CPC cases
#fyaml = 'cice6rest_files_SFS_RTOFS.yaml'
iconc  = 0  # insert iconc from NRT NSIDC
ithkn  = 0  # insert ithkn clim
hsnow  = 1  # insert snow depth clim
snitd  = 1  # snow redistribution over ice (ITDrdg - requires extra tracers in restart) 
sstmom = 0  # adjust MOM6 SST to bring closer to Tfreeze(S) under sea ice, will modify MOM.res.nc 

# enmb option :  For GFS ensembles with different restart files under mem000, mem001, ... dirs
#                Not needed otherwise
# Provide rdate_in (restart date YYYYMMDD) to override YAML value
# provide rhr_in (restart hour, 0,..., 24) to override YAML value

parser = argparse.ArgumentParser()
parser.add_argument("--rdate_in", type=int, default=None,
                    help="Override restart date rdate_in from YAML, YYYYMMDD")
parser.add_argument("--rhr_in", type=int, default=None, 
                    help="Override restart hour rhr_in from YAML, 0,...,24")
parser.add_argument("--iconc", type=int, 
                    choices=[0,1], help=f"insert ice conc NRT NSIDC, default={iconc}")
parser.add_argument("--ithkn", type=int, 
                    choices=[0,1], help=f"insert ice thickn climatology, default={ithkn}")
parser.add_argument("--hsnow", type=int, 
                    choices=[0,1], help=f"insert snow depth climatology, default={hsnow}")
parser.add_argument("--snitd", type=int, 
                    choices=[0,1], help=f"use snow redistribution over ice, default={snitd}")
parser.add_argument("--regn", help=f"region where restart is being updated, default=global", 
                    choices=['south','north','global'], default="global", type=str)
parser.add_argument("--sst", help=f"adjust SST in MOM restart under sea ice, default={sstmom}",
                    choices=[0,1], default=sstmom, type=int) 
parser.add_argument("--fyaml", 
    help=f"YAML with local directories, filenames, restart dates, default={default_yaml}, inside ./config",
    default=default_yaml,
    type=str)
parser.add_argument("--enmb",
      help=f"For GFS ensembles: number for ensamble run, requires YAML for enmb option to construct paths", 
      default=None, type=int)
args = parser.parse_args()

regn  = args.regn
iconc = args.iconc  if args.iconc  is not None else iconc
ithkn = args.ithkn  if args.ithkn  is not None else ithkn
hsnow = args.hsnow  if args.hsnow  is not None else hsnow
snitd = args.snitd  if args.snitd  is not None else snitd
fyaml = CONFIG_DIR / args.fyaml 
enmb  = args.enmb 
sstmom = args.sst
rdate_in = args.rdate_in
rhr_in   = args.rhr_in

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

# Output Restart dates if missing - same as input
# input dates are deduced from restart input file
# or provided in YAML
def to_int_or_none(val):
  return None if val is None else int(val)

rdate_out = to_int_or_none(config_rest["restart_time"]["rdate_out"])
rhr_out   = to_int_or_none(config_rest["restart_time"]["rhr_out"])
# If command-line argument not provided - read date / hr restart input from YAML
if rdate_in is None:
  rdate_in  = to_int_or_none(config_rest["restart_time"]["rdate_in"])
if rhr_in is None:
  rhr_in    = to_int_or_none(config_rest["restart_time"]["rhr_in"])

# Input restart:
# Where original CICE6 restart file is located:
# first, check if there is time offset in sea ice rest. file:
dlt_hrice = config_rest["restart_time"]["dlt_hrice"]
if dlt_hrice is not None:
  date_hrs  = addhours_date(rdate_in, dlt_hrice)
  rdate_ice = int(date_hrs[:8])   # YYYYMMDD
  hr_ice    = int(date_hrs[8:10]) # hour
else:
  rdate_ice = rdate_in
  hr_ice    = rhr_in

if enmb is None:
  #pthrst_in = config_rest["cice_paths"]["pth_in"].format(rdate_ice=rdate_ice, hr_ice=hr_ice)
  pthrst_in = config_rest["cice_paths"]["pth_in"].format(rdate_ice=rdate_ice, hr_ice=hr_ice, IC_DIR=os.environ.get("IC_DIR"))
  flrst_in  = config_rest["rest_names"]["flrst_in"]
else:
  pthrst_in = config_rest["cice_mem_paths"]["pth_in"].format(enmb=f"{enmb:03d}")
  flrst_in  = config_rest["rest_names"]["mem"]["flrst_in"]

# Create flrst_in if it is not the actual name in YAML:
if rdate_in is not None:
  flrst_in = change_base_template(flrst_in, rdate_in, rhr_in, None, err_break = False)

if flrst_in is None:
  raise RuntimeError("Input restart file name is missing in YAML")
else:
  flrst_in_start  = flrst_in

# Derive restart time from file if not specified in YAML
if rdate_in is None or rhr_in is None:
  with xarray.open_dataset(os.path.join(pthrst_in, flrst_in)) as ds:
    try:
      year  = ds.attrs["myear"]
      month = ds.attrs["mmonth"]
      day   = ds.attrs["mday"]
      sec   = ds.attrs["msec"]
    except KeyError as err:
      raise RuntimeError(f"Missing expected date attribute in restart file: {err}")

    rdate_in = int(year*10000 + month*100 + day)
    rhr_in   = sec // 3600

# By now rdate_in should have been identified:
assert rdate_in is not None, f"rdate_in was not identified, check {fyaml} or {flrst_in}"

if rdate_out is None:
  rdate_out = rdate_in
if rhr_out is None:
  rhr_out = rhr_in


print(f"Restart date / time input   :   {rdate_in} : {rhr_in:02d}")
print(f"Restart date / time output  :   {rdate_out} : {rhr_out:02d}")


# Output restart:
#pthrst_out = config_rest["cice_paths"]["pth_out"].format(regn=regn)
#flrst_out  = config_rest["rest_names"]["flrst_out"]
if enmb is None:
  pthrst_out = config_rest["cice_paths"]["pth_out"]
  flrst_out = config_rest["rest_names"]["flrst_out"]
else:
  pthrst_out = config_rest["cice_mem_paths"]["pth_out"].format(enmb=f"{enmb:03d}")
  flrst_out = config_rest["rest_names"]["mem"]["flrst_out"]

# If restart out is not specified, check if template has been provided 
# to use in constructing file name:
if enmb is None:
  flrst_out_tmp = config_rest["rest_names"]["flrst_tmp"]
else:
  flrst_out_tmp = config_rest["rest_names"]["mem"]["flrst_tmp"]

flrst_out_start = flrst_out

# Temporary file names passed along the processes with suffixes being added
# after inserting fields .iconc or .hsnow etc.
# construct temporary file name 
# from restart input if provided otherwise use default file name
# use template file name if provided:
if flrst_out_tmp is not None:
  fltmp_base = change_base_template(flrst_out_tmp, rdate_out, rhr_out, flrst_in)
else:
  # Template not provided, construct from input restart but replace date:
  if flrst_in_start is not None:
    flrst_time_new = flname_replace_date(flrst_in_start, rdate_out, rhr_out)
    fltmp_base = os.path.splitext(flrst_time_new)[0]
  else: 
    fltmp_base = f'cice_restart.{rdate_out}.{rhr_out:02d}'  # default

print(f"fltmp_base={fltmp_base}")

os.makedirs(pthrst_out, exist_ok=True)

print(f"Restart input  directory:\n  {pthrst_in}")
print(f"Restart input  file:\n  {flrst_in_start}")
print(f"Restart output directory:\n  {pthrst_out}")
print(f"Restart output file (if None, will be constructed using {fltmp_base}):\n  {flrst_out_start}")
print(" ==== START INSERTION ====\n\n")

# Insert ice concentration and ice thickness if defined:
if ithkn == 1 or iconc == 1:
  # Create temporary restart if other changes are following:
  # Use it as the final name if flrst_out is not provided
  flrst_out = f'{fltmp_base}.iconc.nc'
  if ithkn == 1:
    flrst_out = f'{fltmp_base}.iconc_ithkn.nc'

  if hsnow == 0 and snitd == 0:
    # No other restart updates, final output name:
    if flrst_out_start is not None:
      flrst_out = flrst_out_start 
      
  print(f"Driver: running ithkn iconc --> CICE6 restart")
  print(f"{flrst_in} --> {flrst_out}")

  script_run = CODE_DIR / "insert_iconc_ithkn_cice6rest_global.py"
  cmd_iconc = [
    "python", str(script_run),
    "--rdate", str(rdate_in),
    "--rhr", str(rhr_in),
    "--rdate_out", str(rdate_out),
    "--rhr_out", str(rhr_out),
    "--pth_in", str(pthrst_in),
    "--flrst_in", str(flrst_in),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
    "--ithkn", str(ithkn),
    "--regn", str(regn),
    "--fyaml", str(fyaml),
  ]

  subprocess.run(cmd_iconc, check=True)

# Insert snow depth:
if hsnow == 1:
  # Update input/output restart files:
  # Output directory keep the same
  if ithkn == 1 or iconc == 1:
    # Continued, use previous restart names:
    flrst_in   = flrst_out
    pthrst_in  = pthrst_out 
    fltmp_base = os.path.splitext(flrst_out)[0]
  else:
    flrst_in = flrst_in_start

  flrst_out = f'{fltmp_base}.hsnow.nc'
  if snitd == 0:
    # Final update:
    if flrst_out_start is not None:
      flrst_out = flrst_out_start
   
  print(f"Driver: inserting hsnow --> CICE6 restart")
  print(f"{flrst_in} --> {flrst_out}")

  script_run = CODE_DIR / "insert_hsnow_cice6rest_global.py"
  cmd_hsnow = [
    "python", str(script_run),
    "--rdate", str(rdate_in),
    "--rhr", str(rhr_in),
    "--rdate_out", str(rdate_out),
    "--rhr_out", str(rhr_out),
    "--pth_in", str(pthrst_in),
    "--flrst_in", str(flrst_in),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
    "--regn", str(regn),
    "--fyaml", str(fyaml),
  ]
    
  subprocess.run(cmd_hsnow, check=True)

if snitd == 1:
  # Need to turn on flags in ice_in to use ITS and snow physics in CICE6
  if ithkn == 1 or iconc == 1 or hsnow == 1:
    # Continued:
    flrst_in = flrst_out
    pthrst_in = pthrst_out
    fltmp_base = os.path.splitext(flrst_out)[0]
  else:
    flrst_in = flrst_in_start
    
  flrst_out = f"{fltmp_base}.snphys.nc"
  # Final update:
  if flrst_out_start is not None:
    flrst_out = flrst_out_start

  print(f"Driver: adding snow ITD and physics variables --> CICE6 restart")
  print(f"{flrst_in} --> {flrst_out}")

  script_run = CODE_DIR / "add_snowvar_snowphys_cice6rest_global.py"
  cmd_snitd = [
    "python", str(script_run),
    "--rdate", str(rdate_in),
    "--rhr", str(rhr_in),
    "--rdate_out", str(rdate_out),
    "--rhr_out", str(rhr_out),
    "--pth_in", str(pthrst_in),
    "--flrst_in", str(flrst_in),
    "--pth_out", str(pthrst_out),
    "--flrst_out", str(flrst_out),
    "--flrst_out", str(flrst_out),
  ]

  subprocess.run(cmd_snitd, check=True)


if sstmom == 1: 
  # Adjust MOM6 SST  under sea ice
  nsec = int(rhr_in * 3600)
  # Can provide MOM6 restart name in fyaml 
  #flmom_in = f"{rdate_in}.{nsec:06d}.MOM.res.nc"
  flmom_in = None
  print(f"Running MOM6 SST update, YAML file = {fyaml}, MOM orig restart = {flmom_in}")

  script_run = CODE_DIR / "correct_surfT_mom6_restart.py"

  assert fyaml is not None, "FYAML name has to be provided"

  cmd_sst = [
    "python", str(script_run),
    "--fyaml", str(fyaml)
  ]

  if flmom_in is not None:
    cmd_sst.extend([
      "--flmom_in", str(flmom_in)
  ])

  if flrst_out is not None:
    cmd_sst.extend([
        "--flice", str(flrst_out)
  ])

  subprocess.run(cmd_sst, check=True)


