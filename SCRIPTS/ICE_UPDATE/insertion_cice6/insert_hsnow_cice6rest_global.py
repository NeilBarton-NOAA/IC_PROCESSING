"""
  Modify snow depth on sea ice in the CICE6 restart file
  by direct insertion of snow depth fields

  Use for modifying Arctic and / or Antarctic 

  Assumed that all non-nan non-zero values are inserted 
  to the grid values where aice > 0

  Snow is distributed across the thikn. categories proportional 
  to the aice (ice partial area)

  Input: snow thickness over sea ice (aggregated over thickness cats)
         on UFS grid

  author: Dmitry Dukhovskoy, NOAA NWS OMD 2026
"""
import os
import numpy as np
import sys
import xarray
import yaml
import argparse
from pathlib import Path
import sys
                      
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from MyPython.mod_time import datevec, rdate2datenum
from MyPython.mod_cice6_utils import (
    ice_enthalpy_BL99, 
    ice_enthalpy_to_temp, 
    adjust_snow_freeboard,
    snow_ice_freeboard, 
    adjust_ice_freeboard
)

CONFIG_DIR = ROOT / "config"
default_yaml = next(CONFIG_DIR.glob("*.yaml")).name


def extract_suffix(fname):
  parts = fname.split('.')
  # must end with .nc, so suffix is the part before that
  if len(parts) > 2 and parts[-1] == 'nc':
    suffix = parts[-2]
    if not suffix.isdigit():
      return suffix
  return None

def check_print(aice, aicen_new, vicen_new, vsnon_new, qsnon, qsnon_new, ncat):
  # Check hice(n) as it is caclulated in icepack_therm_vertical.F90
  # hice(n) = vice(n) / aice(n) 
  print(' =========  ICE  =========')
  for k in range(1,ncat+1):
    aice_n = aicen_new[k-1,:].squeeze()
    vice_n = vicen_new[k-1,:].squeeze()
    hice_n = np.divide(vice_n, aice_n, out=np.zeros_like(aice), where=aice_n != 0)
    jmin, imin = np.unravel_index(hice_n.argmin(), hice_n.shape)
    jmax, imax = np.unravel_index(hice_n.argmax(), hice_n.shape)
    print(f"Cat {k}, j={jmin}, i={imin}, min hice(n): {np.nanmin(hice_n)}, "+\
          f"aice(n): {aice_n[jmin,imin]}, vice(n): {vice_n[jmin,imin]}")
    print(f"  j={jmax}, i={imax}, max hice(n): {np.nanmax(hice_n)}, "+\
          f"aice(n): {aice_n[jmax,imax]}, vice(n): {vice_n[jmax,imax]}")

  print(' =========  SNOW =========')
  for k in range(1,ncat+1):
    aice_n = aicen_new[k-1,:].squeeze()
    vsno_n = vsnon_new[k-1,:].squeeze()
    hsno_n = np.divide(vsno_n, aice_n, out=np.zeros_like(aice), where=aice_n != 0)
    jmin, imin = np.unravel_index(hsno_n.argmin(), hsno_n.shape)
    jmax, imax = np.unravel_index(hsno_n.argmax(), hsno_n.shape)
    print(f"Cat {k}, j={jmin}, i={imin}, min hsnow(n): {np.nanmin(hsno_n)}, "+\
          f"aice(n): {aice_n[jmin,imin]}, vsno(n): {vsno_n[jmin,imin]}")
    print(f"  j={jmax}, i={imax}, max hsnow(n): {np.nanmax(hsno_n)}, "+\
          f"aice(n): {aice_n[jmax,imax]}, vsno(n): {vsno_n[jmax,imax]}")

  print(' ======== snow enthalpy =======')
  for k in range(1,ncat+1):
    qtot_old = np.nansum(qsnon[k-1,:])
    qtot_new = np.nansum(qsnon_new[k-1,:])
    print(f"Cat {k}, old snow enth={qtot_old:.4e} new snow enth={qtot_new:.4e}")

  print(" ")

  return

def insert_hsnow(regn_wrk, ds_out, config_rest, dnmbR, dnmbN, pthrest, flrst_in):
  """
    regn_wrk    - region name in this loop: north, south 
    ds_out      - xarray data_set: empty or not if cycled over 2 regions
    config_rest - YAML object with local directories, files, restart fields, etc.
    dnmbR       - input restart date datenumb format
    dnmbN       - output restart date datenumb format
    pthrest     - restart directory
    flrst_in    - input restart file to be modified    
  """
  # CICE parameters:
  puny      = 1.e-11
  c0        = 0.0
  c1        = 1.0
  c2        = 2.0
  p5        = 0.5
  Lsub      = 2.835e6    # latent heat sublimation fw (J/kg)
  Lvap      = 2.501e6    # latent heat vaporization fw (J/kg)
  Lfresh    = Lsub - Lvap # latent heat of melting of fresh ice (J/kg)
  cp_ice    = 2106.       # specific heat of fresh ice (J/ kg/K)
  rhos      = 330.        # density of snow (kg/m3)
  hs_min    = 1.e-4       # min snow thickness for computing Tsno (m)
  nsal      = 0.407
  msal      = 0.573
  min_salin = 0.1      # threshold for brine pocket treatment
  saltmax   = 3.2        # max S at ice base
  hg        = 1.e20    # bad values, land mask, etc.
  nslyr     = 1   # snow layers
  Tmin      = -100.   # minimum snow T
  # Check ice_in what ITD is used
  # 0.00, 0.64, 1.39, 2.47, 4.57
  hicat = np.array([0., 0.64, 1.39, 2.47, 4.57, 50.])

  # in CICE6, snow surf max T = 0C
  # Set surface snow T < 0 to prevent rapid snow melt during the first time steps
  # This mainly applies for summer months
  # Tsfc = Tsnow in the 1 layer --> change qsnon(1)
  Tsfc_max = -1.0   
  rho_ice  = 917. 
  hsnow_max = 500.   # to avoid very thick hsnow / ice_area which will cause picard iteration crush
  if regn_wrk == 'south':
    rho_ocean = 1025.
  elif regn_wrk == 'north':
    rho_ocean = 1025.

  # Checks:
  assert nslyr == 1, f"Code needs to be modified for nslyr>1, nslyr={nslyr}"
  assert Tsfc_max <= 0., f"Tsfc_max={Tsfc_max} has to be <=0"

  # Restart time: in
  yrR,mmR,ddR,hrR = datevec(dnmbR, round_hrs=True)[:4]
  nsecR = hrR*3600
    
  # Restart time: out
  yrN, mmN, ddN, hrN = datevec(dnmbN, round_hrs=True)[:4]
  nsecN = hrN*3600
  
  # Snow depth climatology, Interpolated fields mesh025:
  varnm = 'snow_depth'
  if regn_wrk == 'south':
    pthsnow = config_rest["target_paths"]["hsnow_south"]["path"]
    flhsn   = config_rest["target_paths"]["hsnow_south"]["file"]
  elif regn_wrk == 'north':
    pthsnow = config_rest["target_paths"]["hsnow_north"]["path"]
    flhsn   = config_rest["target_paths"]["hsnow_north"]["file"]

  dflhsn = os.path.join(pthsnow,flhsn)

  print(f"Reading interpolated hsnow {dflhsn}")
  with xarray.open_dataset(dflhsn) as ds_snow:
    HSi = ds_snow[varnm].isel(time=mmN-1).data.squeeze()
    LON = ds_snow['lon'].data
    LAT = ds_snow['lat'].data
    units = ds_snow[varnm].attrs.get('units', None)
    if units is not None:
      print(f"'snow_depth' units: {units}")
      hunits = units
    else:
      print("No 'units' attribute found for 'snow_depth', use default: {hunits}")

  units_m = hunits == 'm'
 
  if ds_out is None: 
    dflrst_in = os.path.join(pthrest, flrst_in)
    print(f"Reading restart: {dflrst_in}")
    ds_in = xarray.open_dataset(dflrst_in)
    ds_out = ds_in.copy(deep=True)
    ds_in.close()
  else:
    print(f"Initialization of ds_out is skipped, continue with existing ds_out ...")

  # Input values:
  aicen  = ds_out['aicen'].data   # partial area by cats
  vsnon  = ds_out['vsnon'].data   # snow vol per m2 of ice area
  qsnon  = ds_out['qsno001'].data # snow enthalpy by cats for 1 snow layer
  vicen  = ds_out['vicen'].data   # ice vol per unit area of grid cell m3/m2
  tsfcn  = ds_out['Tsfcn'].data   # snow/ice surface T
  qicen1 = ds_out['qice001'].data # ice enthalpy, lr 1 surface
  sicen1 = ds_out['sice001'].data # ice S, layer 1
  apndn  = ds_out['apnd'].data    # the fraction of the pond of ice area, for each cat
  hpndn  = ds_out['hpnd'].data    # depth of the ponds in a cell, by cats
  ncat, jdim, idim = vsnon.shape

  # Aggregated ice partial area:
  aice = np.sum(aicen, axis=0).squeeze()

  # Select points to insert:
  Jins, Iins = np.where((HSi > puny) & (~np.isnan(HSi)) & (aice > puny))
  Xins = LON[Jins,Iins]
  Yins = LAT[Jins,Iins]
  npnts = len(Jins)

  print(f"Found {npnts} points for insertion, min/max lat={np.min(Yins):.1f}/{np.max(Yins):.1f}"
         f" lon={np.min(Xins):.1f}/{np.max(Xins):.1f}")

  # Note qsnon, qice < 0 !
  aicen_new  = aicen.astype(ds_out['aicen'].dtype).copy()
  vsnon_new  = vsnon.astype(ds_out['vsnon'].dtype).copy()
  qsnon_new  = qsnon.astype(ds_out['qsno001'].dtype).copy()
  qicen1_new = qicen1.astype(ds_out['qice001'].dtype).copy()
  apndn_new  = apndn.astype(ds_out['apnd'].dtype).copy()
  hpndn_new  = hpndn.astype(ds_out['hpnd'].dtype).copy()
  tsfcn_new  = tsfcn.astype(ds_out['Tsfcn'].dtype).copy()
  vicen_new  = vicen.astype(ds_out['vicen'].dtype).copy()

  dvol_sum = 0.
  print("Snow depth insertion ...")
  for ipp in range(npnts):
    if ipp%10000 == 0:
      print(f"   {ipp/npnts*100.:.2f}% done ...")
    j0 = Jins[ipp]
    i0 = Iins[ipp]

    ai  = aice[j0,i0]           # aggregated ice partial area 
    ain = aicen[:,j0,i0]        # partial areas by cats
    vin = vicen[:,j0,i0]        # ice volume per unit grid cell area by cats
    vsn = vsnon[:,j0,i0]        # snow volume per unit grid-cell area m2
    apn = apndn[:,j0,i0]        # pond fractional area of ice (level) area
    hpn = hpndn[:,j0,i0]        # pond depth
    tsn = tsfcn[:,j0,i0]        # snow/ice surface T by cats

    # New ice snow thickness over sea ice:
    if units_m:
      hsn_new = HSi[j0,i0]        # m of snow over sea ice
    else:
      hsn_new = HSi[j0,i0]*0.01   # m of snow over sea ice 

    #check_pnt =  LAT[j0,i0] > 70 and ai > 1
    #if check_pnt:
    #  print(f"checking: lon={LON[j0,i0]:.2f}, lat={LAT[j0,i0]:.2f}")
    #  print(f"i0={i0} j0={j0}, ai={ai:.3f}, hsn_new={hsn_new:.5f}, HSi={HSi[j0,i0]:.5f}")    

    # This should not happen (target HSi assumed to be m3/m2_ice) just to make sure
    # there is no very thick snow --> possible crash in picard iteration for very thick snow or ice
    if hsn_new > hsnow_max:
      hsn_new = hsnow_max

    # ice conc should not change except for a few locaitons to adjust snow load across cats:
    ain_new = ain.copy()
    
    # Distribute new snow depth evenly by cats in snow vol m3/m2:
    vstot_new = hsn_new * ai
    if hsn_new <= hs_min or ai < puny:
      vsn_new = np.zeros_like(ain)
    else:
      # Distribute across cats proportionally to iconc:
      wts = ain/ai
      vsn_new = vstot_new * wts

    # Update snow enthalpy: J/m3  
    # snow enthalpy should be: qsn_min <= qsn <= qsn_max
    # In theory, qsn_max = -rhos*Lfresh (latent heat of metling at 0C)
    # Make it a little lower to keep snow from melting right away
    # In general, snow enth. = enth(Tsfcn) if Tsfcn <=0
    qsn = qsnon[:,j0,i0]                                      # enthalpy, J/m3 < 0
    qsn_min = -rhos * Lfresh + (Tmin + 0.01) * cp_ice * rhos  # enth. of the coldest possible snow
    qsn_max = -rhos * Lfresh + Tsfc_max * cp_ice * rhos       # Tsfc_max <= 0
    qsn_tsf = -rhos * Lfresh + tsn * cp_ice * rhos            # enth. for surf temp
    qT0 = -Lfresh*rhos                                        # enth. of pure snow at 0C

    # Update enthalpy of snow and enforce physical constraints
    # Limit qsn to [qsn_min, qsn_max]
    qsn_new = np.clip(qsn, qsn_min, qsn_max)

    # It should not exceed enthalpy implied by surface temperature (qsn </= qsn_tsf)
    # This is true for 1 snow layer
    qsn_new = np.minimum(qsn_new, qsn_tsf)

    # No ice --> no snow enthalpy
    qsn_new = np.where(ain <= puny, 0., qsn_new)

    # Zero snow volume --> zero enthalpy
    qsn_new = np.where(vsn_new <= 0., 0., qsn_new)

    # Tsfcn should match snow enthalpy in layer 1
    # Update Tsfcn if needed:
    # for dry snow (T<=0):
    tsn_new = (qsn_new + rhos*Lfresh) / (cp_ice*rhos)
    tsn_new = np.where(abs(qsn_new) < puny, 0., tsn_new)
    for ik in range(ncat):
      if abs(qsn_new[ik]) < puny:
        continue
      assert tsn_new[ik] <= 0., f"check qsn_new={qsn_new[ik]:.4e} --> tsn_new={tsn_new[ik]}"

    # Update ice enthalpy in the surface layer to prevent rapid snow melt
    # if qice > qsnow, this is particularly important for 
    # no snow --> snow cases during summer
    sin = sicen1[:,j0,i0]
    qin = qicen1[:,j0,i0] 
    qin_new = ice_enthalpy_BL99(tsn_new, sin)
    qin_new = np.minimum(qin_new, qin)
    # Check ice:
    Tice_new = ice_enthalpy_to_temp(qin_new,sin)
    mu_ice = 0.054  # liquidus ratio btw frz T and salinity of brine, [deg/ppt] BL99
    Tice_melt = -mu_ice * sin
    if np.any(Tice_new >= Tice_melt):
      print(f"j0={j0}, i0={i0}, ice T exceeds melting T")
      for kcat in range(len(Tice_new)):
        print(f"cat={kcat+1}: Tice_new={Tice_new[kcat]:.4f}  Tmelt={Tice_melt[kcat]:.4f}")
      raise Exception("ERR Updating ice enthalpy layer 1")

    vtot_init = np.nansum(vsn)
    vtot_new  = np.nansum(vsn_new)
    #print(f"tot vsnon change = {vtot_new-vtot_init}") 

    # Pond fractional area and depth
    # For now, make it 0 to prevent rapid snow melting when apnd >> 0
    # Note if apnd is changed and > 0, need to adjust hpnd to maintain nonnegative freeboard
    # see: icepack_meltpond_lvl.F90 as an example
    apn_new = apn * 0.
    hpn_new = hpn * 0.

    # Snow-ice interface should be at or above the sea level
    # if it is below, snow will be melted instantly to bring the interface to sea level
    # First, if needed - try to redistribute excess snow across ice thikn. cats:
    vsn_new = adjust_snow_freeboard(vin, vsn_new, ain, rho_ice=rho_ice, \
                              rho_snow=rhos, rho_ocean=rho_ocean)

    # Check if the snow-ice interface in all cats is above the sea level
    # if not, try to slightly modify ice in the cat where needed
    # Note this will slightly change ice volume (ice thickness)
    # skip this step if ice vol has to be preserved
    ice_frb = snow_ice_freeboard(vin, vsn_new, ain, \
                      rho_ice=rho_ice, rho_snow=rhos, rho_ocean=rho_ocean)  
    if np.any(ice_frb < 0.):
      vin_new, ain_new, vsn_new = adjust_ice_freeboard(vin, vsn_new, ain, hicat, \
                      rho_ice=rho_ice, rho_snow=rhos, rho_ocean=rho_ocean, fdebug=False) 
    else:
      vin_new = vin.copy()

    #if check_pnt:
    #  print(f"Check: vsn_new={vsn_new}")
    #  print(f"       sum vsn_new = {np.sum(vsn_new)}")

    # Update:
    dvol_sum = dvol_sum + (vtot_new-vtot_init)
    aicen_new[:,j0,i0]  = ain_new
    vsnon_new[:,j0,i0]  = vsn_new
    qsnon_new[:,j0,i0]  = qsn_new
    qicen1_new[:,j0,i0] = qin_new
    apndn_new[:,j0,i0]  = apn_new
    hpndn_new[:,j0,i0]  = hpn_new
    tsfcn_new[:,j0,i0]  = tsn_new
    vicen_new[:,j0,i0]  = vin_new

    diff = np.nansum(vsn_new - vsnon[:, j0, i0])
    diff2 = np.nansum(vsn_new -vsn)
    diff3 = np.nansum(vsnon[:,j0,i0] - vsnon_new[:,j0,i0])
    if diff == 0 and abs(diff2) > 0:
      print(f"No change at {j0},{i0}, expected diff={diff2}")

    if diff3 == 0 and abs(diff2) > 0:
      print(f"No change in the arrays at {j0},{i0}, expected diff={diff2}")

  # Checking:
  print(f"Snow vol change: dvol_sum = {dvol_sum}")
  total_vsnon_init = np.nansum(vsnon)
  total_vsnon_new  = np.nansum(vsnon_new)
  print(f"Tot snow volume change (m3/m2): {total_vsnon_new - total_vsnon_init}")

  # Update data set:
  #ds_out = ds_out.assign(vsnon=vsnon_new, qsno001=qsnon_new)
  ds_out['aicen'].values[:]   = aicen_new
  ds_out['vsnon'].values[:]   = vsnon_new
  ds_out['qsno001'].values[:] = qsnon_new
  ds_out['qice001'].values[:] = qicen1_new
  ds_out['apnd'].values[:]    = apndn_new
  ds_out['hpnd'].values[:]    = hpndn_new
  ds_out['Tsfcn'].values[:]   = tsfcn_new
  ds_out['vicen'].values[:]   = vicen_new

  # Sanity checking:
  assert "vsnon" in ds_out and "qsno001" in ds_out, "Missing updated snow fields vsnon and qsno001"
  assert ds_out["vsnon"].shape == vsnon_new.shape, "Check shape of vsnon "
  assert ds_out["qsno001"].shape == qsnon_new.shape, "Check shape of qsnon "

  check_print(aice, aicen_new, vicen_new, vsnon_new, qsnon, qsnon_new, ncat)

  return ds_out


def main():
  rest_hr   = 0
  hunits    = 'cm'

  yrR = mmR = ddR = hrR = None
  yrN = mmN = ddN = hrN = None

  parser = argparse.ArgumentParser()
  parser.add_argument("--rdate", help=f"restart date input file", required=True, type=int)
  parser.add_argument("--rhr", help=f"input file, restart hour = 0, ..., 23, default={rest_hr}", type=int)
  parser.add_argument("--rdate_out", help="output file, restart date", required=True, type=int)
  parser.add_argument("--rhr_out", help="output file, restart hour", required=True, type=int)
  parser.add_argument("--pth_in", help="input restart directory with original file", type=str, required=True)
  parser.add_argument("--flrst_in", help="rest file in", required=True, type=str)
  parser.add_argument("--pth_out", help="output restart directory where new file be dumped", type=str, required=True)
  parser.add_argument("--flrst_out", help="new rest file", required=True, type=str)
  parser.add_argument("--regn", help=f"where hsnow inserted", type=str, required=True, 
                       choices=['north','south','global'])
  parser.add_argument("--fyaml",
    help=f"YAML with local directories, filenames, restart dates, default={default_yaml}, inside ./config",
    default=default_yaml,
    type=str)
  args = parser.parse_args()

  flrst_in      = args.flrst_in 
  flrst_out     = args.flrst_out
  regn          = args.regn   
  rest_date     = args.rdate
  rest_hr       = args.rhr 
  rest_date_out = args.rdate_out 
  rest_hr_out   = args.rhr_out 
  pth_in        = args.pth_in
  pth_out       = args.pth_out
  fyaml         = CONFIG_DIR / args.fyaml
  
  print(f"Restart date input:  {rest_date}:{rest_hr}")
  print(f"Restart date output: {rest_date_out}:{rest_hr_out}")

  # Get date numbers:
  # Input restart file
  dnmbR = rdate2datenum(rest_date*100+rest_hr)  # restart day nmb
  if yrR is None:
    yrR, mmR, ddR, hrR = datevec(dnmbR, round_hrs=True)[:4]
  nsecR = hrR*3600

  # Dates of the output fields in the new restart:
  dnmbN = rdate2datenum(rest_date_out*100+rest_hr_out)
  if yrN is None:
    yrN, mmN, ddN, hrN = datevec(dnmbN, round_hrs=True)[:4]
  nsecN = hrN*3600
   
  with open(fyaml) as ff:
    config_rest = yaml.safe_load(ff)

  def expand_env_vars(d):
      for k, v in d.items():
          if isinstance(v, dict):
              expand_env_vars(v)
          elif isinstance(v, str):
              d[k] = os.path.expandvars(v)

  expand_env_vars(config_rest)

  # Insert hsnow for regions
  ds_out = None
  if regn == 'global':
    for regn_tmp in (['north','south']):
      print(f"  hsnow Processing region: {regn_tmp}")
      ds_out = insert_hsnow(regn_tmp, ds_out, config_rest, dnmbR, dnmbN, pth_in, flrst_in)
  else:
    print(f"  hsnow Processing region: {regn}")
    ds_out = insert_hsnow(regn, ds_out, config_rest, dnmbR, dnmbN, pth_in, flrst_in)

  # Attributes:
  from datetime import datetime
  istep1_val = ds_out.attrs.get('istep1', None)
  ds_out.attrs.update({
      "title": "CICE6 restart with inserted hsnow from climatology gridded fields",
      "source": "insert_hsnow_cice6rest_global.py",
      "istep1": np.int32(istep1_val) if istep1_val is not None else np.int32(0), 
      "myear": np.int32(yrN),
      "mmonth": np.int32(mmN),
      "mday": np.int32(ddN),
      "msec": np.int32(nsecN),
      "history2": f"Modified {flrst_in} {datetime.now().isoformat()}",
      "region": regn,
  })

  if flrst_out is None: 
    if ins_thkn:
      flrst_out = f"cice_restart.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.hsnow.nc"
    else:
      flrst_out = f"cice_restart.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.hsnow.nc"

  dflrst_out = os.path.join(pth_out,flrst_out)
  print(f"Saving CICE restart --> {dflrst_out}")
  ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
  ds_out.close()

if __name__ == "__main__":
  main()
