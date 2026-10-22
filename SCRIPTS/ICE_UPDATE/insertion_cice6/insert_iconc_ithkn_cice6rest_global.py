"""
  Insert ice concentration (sea ice partial area) and 
  (optionally) ice thickness 
  into CICE6 aice fields 
  Using NSIDC interpoalted fields for iconc
  and CryoSat for ice thickness

  The insertion is performed for both Arctic and Antarctic or only one of these regions
  older version: insert_iconc_ithkn_cice6_restart.py performs this for 1 region only

  if ice thickness is opted out, this should be idential to
  insert_iconc_cice6_restart.py

  See python/prepare_cice6/interp_NSIDC_iconc_mesh025.py
                    interp_CryoSat_ithkn_antarct_mesh025.py

  author: Dmitry Dukhovskoy, NOAA NWS OMD 2026
"""
import os
import numpy as np
import xarray
import yaml
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from MyPython.mod_time import datevec, rdate2datenum
from MyPython.mod_mom6 import read_mom6grid
from MyPython.mod_cice6_utils import (
    change_base_template,
    flname_replace_date,
    adjust_thkncats_aice, 
    check_ithkn_cats, 
    sice_lr_cice4,
    ice_enthalpy_BL99
)
 
CONFIG_DIR = ROOT / "config"
default_yaml = next(CONFIG_DIR.glob("*.yaml")).name

def find_adj_icepnts(dlti, aice, i0, j0, puny):
  jdm, idm = aice.shape
  i1 = max(i0 - dlti, 0)
  i2 = min(i0 + dlti, idm - 1)
  j1 = max(j0 - dlti, 0)
  j2 = min(j0 + dlti, jdm - 1)

  A = aice[j1:j2+1, i1:i2+1]
  jmm, imm = np.where(A > puny)

  # No grid points with ice found:
  if jmm.size == 0:
    return np.array([], dtype=int), np.array([], dtype=int)

  jice = jmm + j1
  iice = imm + i1

  return iice, jice

def find_adj_ocnpnts(dlti, aice, i0, j0, puny):
  jdm, idm = aice.shape
  i1 = max(i0 - dlti, 0)
  i2 = min(i0 + dlti, idm - 1)
  j1 = max(j0 - dlti, 0)
  j2 = min(j0 + dlti, jdm - 1)

  A = aice[j1:j2+1, i1:i2+1]
  jmm, imm = np.where(A <= puny)

  # No grid points with no ice found:
  if jmm.size == 0:
    return np.array([], dtype=int), np.array([], dtype=int)

  jocn = jmm + j1
  iocn = imm + i1

  return iocn, jocn

def check_print(aice, aicen_new, vicen_new, vsnon_new, ncat):
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
    print(f"         j={jmax}, i={imax}, max hice(n): {np.nanmax(hice_n)}, "+\
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
    print(f"         j={jmax}, i={imax}, max hsnow(n): {np.nanmax(hsno_n)}, "+\
          f"aice(n): {aice_n[jmax,imax]}, vsno(n): {vsno_n[jmax,imax]}")


  return

def find_varnm(dflithkn, var_opt):
  with xarray.open_dataset(dflithkn) as ds_ithkn:
    for varnm in var_opt:
      if varnm in ds_ithkn.data_vars:
        #print(f"Using variable {varnm}")
        return varnm
        
  raise KeyError("No ice thickness variable name found, check file")

  return

def yaml_path(cfg, **fmt):
  path = cfg["path"].format(**fmt)
  file = cfg["file"].format(**fmt)
  return os.path.join(path, file)

def insert_iconc_ithkn(regn_wrk, ds_out, config_rest, HH, LAT, LON, dnmbR, dnmbN,\
                       ins_thkn, pthrest, flrst_in):
  """
    regn_wrk - region name in this loop: north, south 
    ds_out   - xarray data_set: empty or not if cycled over 2 regions
    config_rest - YAML object with local directories, files, restart fields, etc.
    HH          - topo array
    LAT, LON    - long/ latit arrays
    dnmbR       - input restart date datenumb format
    dnmbN       - output restart date datenumb format
    ins_thkn    - flag: insert thickness  or not with iconc 
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

  # Restart time: in
  yrR,mmR,ddR,hrR = datevec(dnmbR, round_hrs=True)[:4]
  nsecR = hrR*3600
  rdate_in = int(yrR*10000 + mmR*100 + ddR)

  # Restart time: out
  yrN, mmN, ddN, hrN = datevec(dnmbN, round_hrs=True)[:4]
  nsecN = hrN*3600

  fmt = {
      "YR": yrN,
      "MM": mmN,
      "DD": ddN,
      "rdate": rdate_in,
      "regn": regn_wrk
  }

  # Interpolated NSIDC ice conc:
  RMsk = np.where(HH>=0, 0, 1)
  if regn_wrk == 'south':
    RMsk[LAT > -60.] = 0
  elif regn_wrk == 'north':
    RMsk[LAT < 50.] = 0  

  #pthiconc = config_rest["target_paths"]["ice_conc"]["path"].format(YR=yrN)
  #fliconc  = config_rest["target_paths"]["ice_conc"]["file"].format(YR=yrN, MM=mmN, regn=regn_wrk)
  #dfliconc = os.path.join(pthiconc, fliconc)
  dfliconc = yaml_path(config_rest["target_paths"]["ice_conc"], **fmt)
  print(f'Loading target ice conc {dfliconc}')
  with xarray.open_dataset(dfliconc) as dsint:
    AICEint = dsint['ice_conc'].isel(time=ddN-1).squeeze()

  AICEint = np.where(RMsk == 0, np.nan, AICEint)

  # Read ice thickness data if thickness is inserted:
  if ins_thkn:
    if regn_wrk == 'south':
      #pthithkn = config_rest["target_paths"]["ithkn_south"]["path"]
      #flithkn  = config_rest["target_paths"]["ithkn_south"]["file"]
      dflithkn = yaml_path(config_rest["target_paths"]["ithkn_south"], **fmt)
    elif regn_wrk == 'north':
      #pthithkn = config_rest["target_paths"]["ithkn_north"]["path"]
      #flithkn  = config_rest["target_paths"]["ithkn_north"]["file"]
      dflithkn = yaml_path(config_rest["target_paths"]["ithkn_north"], **fmt)

    #dflithkn = os.path.join(pthithkn,flithkn)

    var_opt = ['ice_thkn', 'ithkn', 'hi', 'ice_thickness']
    ithkn_varnm = find_varnm(dflithkn, var_opt)
    print(f"Reading ice thickn varnm='{ithkn_varnm}' for month {mmN} from {dflithkn}")

    with xarray.open_dataset(dflithkn) as ds_ithkn:
      # Assume:
      # ntimes = 1  -> single field
      # ntimes = 12 -> monthly climatology

      ntimes = ds_ithkn.sizes["time"]
      if ntimes == 1:
        ITHKN = ds_ithkn[ithkn_varnm].values.squeeze()
      elif ntimes == 12:
        ITHKN = ds_ithkn[ithkn_varnm].isel(time=mmN-1).values
      else:
        raise RuntimeError(f"Unexpected number of time levels: {ntimes}")

  else:
    # No ice thickness insertion
    ITHKN = np.full_like(HH, np.nan)

  if flrst_in is None:
    #flrst_in = f"cice_model.res.{yrR}{mmR:02d}{ddR:02d}.{nsecR:06d}.nc"
    raise RuntimeError("insert_iconc_ithkn: input restart filename is missing ...")

  if ds_out is None:
    dflrst_in = os.path.join(pthrest, flrst_in)
    print(f"Reading restart: {dflrst_in}")
    ds_in = xarray.open_dataset(dflrst_in)
    ds_out = ds_in.copy(deep=True)
    ds_in.close()
  else:
    print(f"Input restart is skipped, continue with existing ds_out")


  # Insert aice and distribute/reduce proportionally over ice cats.:
  assert nslyr == 1, f"Code needs to be modified for nslyr>1, nslyr={nslyr}"
  aicen = ds_out['aicen'].data  # partial area by cats
  vicen = ds_out['vicen'].data  # ice vol per m2 of grid cell by cats
  vsnon = ds_out['vsnon'].data  # snow vol per m2 of ice area
  qsnon = ds_out['qsno001'].data  # snow enthalpy by cats for 1 snow layer
  Tsfcn = ds_out['Tsfcn'].data  # surface T in each cat.
  ncat, jdim, idim = vsnon.shape

  sice = {}
  qice = {}
  nilrs = 7   # ice layers
  for i in range(1, nilrs+1):
    varnum = f"{i:03d}" 
    sice[varnum] = ds_out[f"sice{varnum}"].data
    qice[varnum] = ds_out[f"qice{varnum}"].data

  # Aggregated ice partial area:
  aice = np.sum(aicen, axis=0).squeeze()

  # Select points to insert:
  mask = (RMsk > 0) & np.isfinite(AICEint)
  Jins, Iins = np.where(mask)
  Xins = LON[Jins, Iins]
  Yins = LAT[Jins, Iins]
  npnts = Jins.size

  npnts_ithkn = 0
  if ins_thkn:
    mask_ithkn = (RMsk > 0) & np.isfinite(ITHKN)
    npnts_ithkn = np.count_nonzero(mask_ithkn)

  print(f"iconc insertion: {npnts} pnts, min/max lat={np.min(Yins):.1f}/{np.max(Yins):.1f}"
         f" lon={np.min(Xins):.1f}/{np.max(Xins):.1f}")
  if ins_thkn:
    print(f"ithkn insertion: {npnts_ithkn} pnts ")

  # Ice thickness cats
  # for kitd = 1 - linear remapping
  # and kcatbound = 0 , the lower cat. thickness values:
  # 0.00, 0.64, 1.39, 2.47, 4.57
  hicat = np.array([0., 0.64, 1.39, 2.47, 4.57, 50.])
  dhi_min = 0.01  # min diff between cat ice thicknesses from 2 adjacent cats
  hcat_indx = np.arange(1,ncat+1)

  # Note qsnon, qice < 0 !
  vsnon_new = vsnon.astype(ds_out['vsnon'].dtype).copy()
  qsnon_new = qsnon.astype(ds_out['qsno001'].dtype).copy()
  aicen_new = aicen.astype(ds_out['aicen'].dtype).copy()
  vicen_new = vicen.astype(ds_out['vicen'].dtype).copy()
  Tsfcn_new = Tsfcn.astype(ds_out['Tsfcn'].dtype).copy()
  qicen_new = {}
  sicen_new = {}
  for i in range(1, nilrs+1):
    varnum = f"{i:03d}"
    qicen_new[varnum] = qice[varnum].astype(ds_out[f"qice{varnum}"].dtype).copy()
    sicen_new[varnum] = sice[varnum].astype(ds_out[f"sice{varnum}"].dtype).copy()

  Tfrz = -1.86243522  # ocean freez. T
  Tsfc_max = -0.1  # max surf temp
  dvol_sum = 0.
  dlti = 2         # N of +/- i , j indices to search for ice pnts around i0,j0
  vitot_min = 0.05    # min total ice vol m3/m2_grid, when ai_old = 0 --> ai_new > 0
  hitot_min = 0.1    # mean ice thickness over ice area: = sum(hice(n)*aice(n)) / sum(aice(n)) 
  hsnow_min = 0.01   # min snow thickness for noice --> ice case, this is m3/m2_ice 
  hsnow_max = 500.   # to avoid very thick hsnow / ice_area which will cause picard iteration crush

  if ins_thkn:
    print("Ice concentration and ice thickness insertion ...")
  else:
    print("Ice concentration insertion ...")
  for ipp in range(npnts):
    if ipp%10000 == 0:
      print(f"   {ipp/npnts*100.:.2f}% done ...")
    j0 = Jins[ipp]
    i0 = Iins[ipp]

    # Note hsnow = vsn / aice for aice > 0
    # for cat n: vsn(n) = hsnow(n) * aice(n) 
    ai_old  = aice[j0,i0]       # aggreageted ice partial area 
    ain_old = aicen[:,j0,i0]    # partial areas by cats
    vsn_old = vsnon[:,j0,i0]    # snow volume per unit grid-cell area m2
    vin_old = vicen[:,j0,i0]    # ice volume per unit grid-cell area m2 
    tsfcn_old = Tsfcn[:,j0,i0]  # surf T

    # Distribute new iconc proportionally by cats in snow vol m3/m2:
    ai_new = AICEint[j0,i0]
    if ai_new <= puny:
      ai_new = 0.
   
    if ai_new > 1.:
      ai_new = 1.

    # iconc change:
    if ai_new < puny:
      ain_new = ain_old * 0.
    else:
      if ai_old < puny:
        # all new iconc in cat 1:
        ain_new = ain_old * 0.
        ain_new[0] = ai_new
      else:
        cff = ai_new/ai_old
        ain_new = ain_old * cff
        
    # Adjust small truncation errors:
    sum_ain = np.sum(ain_new)
    if (sum_ain > 1.0) and (sum_ain - 1.0 < 1e-12):
      ain_new = ain_new / sum_ain - 1.e-12
    ain_new = np.where(ain_new < puny, 0., ain_new)
   
    assert np.sum(ain_new) <= 1., f"Check ain_new: sum>1: {np.sum(ain_new)}"
    assert np.min(ain_new) >= 0., f"Check ain_new: min val < 0 {np.min(ain_new)}"

    # Update surface T and ice vol / ice thickness
    tsf_new = None
    iice = jice = iocn = jocn = None
    vin_new = None
    aice_case = None
    if ai_old <= puny and ai_new > puny:
      # Case: no ice --> ice, created ice in the grid cell
      aice_case = "noice2ice"

      # Update Tsfcn
      # Find N closest ice points:
      iice, jice = find_adj_icepnts(dlti, aice, i0, j0, puny) 
      if len(iice) == 0:
        # no ice pnt adjacent to j0,i0:
        tsf_new = Tsfcn[:,j0,i0]*0.0 + Tsfc_max
        vin_new = hitot_min * ain_new
      else: 
        # where ice - <= Tmax, where no ice = Tfrz
        Tsf_adj = Tsfcn[:,jice,iice]
        tsf_new  = np.nanmean(Tsf_adj, axis=1)
        tsf_new  = np.where(tsf_new > Tsfc_max, Tsfc_max, tsf_new)
        tsf_new  = np.where(ain_new < puny, Tfrz, tsf_new) 

        Vice_adj = vicen[:,jice,iice]
        vin_new = np.nanmean(Vice_adj, axis=1)
        vitot_new = np.max([np.sum(vin_new), vitot_min]) # total ice vol/grid area m3/m2 =grid mean ice thkn, m
        wt = ain_new / np.sum(ain_new)
        vin_new = vitot_new * wt

    elif ai_old > puny and ai_new > puny:
      # Case: ice --> updated ice conc
      aice_case = "ice2ice"
      tsf_new = np.where(tsfcn_old > Tsfc_max, Tsfc_max, tsfcn_old)
      tsf_new  = np.where(ain_new < puny, Tfrz, tsf_new)

      # Try to preserve mean ice thkn over ice:
      vitot_old = np.sum(vin_old)                          # m3/m2_cell or mean ice thkn over grid cell
      hitot_old = vitot_old / ai_old                       # m3/m2_ice, mean ice thkn over ice
      vin_new   = np.max([hitot_old,hitot_min]) * ain_new  # m3/m2_cell or grid-cell mean ice thickness, m

    elif ai_old > puny and ai_new < puny:
      # Case: ice --> no ice
      aice_case = "ice2noice"
      # Copy tsfc from adj grid cells
      iocn, jocn = find_adj_ocnpnts(dlti, aice, i0,j0, puny)
      if len(iocn) == 0:
        # no ocn pnts:
        tsf_new = Tsfcn[:,j0,i0]*0.0 + Tfrz
      else:
        Tsf_adj = Tsfcn[:,jocn,iocn]
        tsf_new = np.nanmean(Tsf_adj, axis=1)
        tsf_new = np.where(tsf_new > Tsfc_max, Tsfc_max, tsf_new)
            
      vin_new = vin_old*0.0

    elif ai_old < puny and ai_new < puny:
      # Case: no ice --> no ice
      aice_case = "noice2noice"
      tsf_new = Tsfcn[:,j0,i0]
      vin_new = vin_old*0.0
    
    else:
      raise Exception(f"Unexpected case for ai_old={ai_old} and ai_new={ai_new}")  

    # Should not happen but Checking if any NaN occur:
    if np.isnan(tsf_new).any():
      tsf_new = np.where(np.isnan(tsf_new), Tsfc_max, tsf_new)
    if np.isnan(vin_new).any():
      vin_new = np.where(np.isnan(vin_new), 0., vin_new)

    # Check that ice vol is correctly distributed across the ice thickn. cats:
    # If not - distribute across ice cats conserving aice and vice 
    # and matching ice cats
    ain_min = 1.e-8    # lower bound of ain(n) to avoid zeros
    vice_clim = 0.     # ice vol / m2_cell from clim (CryoSat is cell mean ice thickn)
    vice_old = 0.      # ice vol / m2 _cell from old restart
    if ins_thkn:
      vice_clim = ITHKN[j0,i0] 
      if np.isnan(vice_clim):
        vice_clim = 0.

    vice_old = np.sum(vin_old)
    vice_new = np.sum(vin_new)

    if vice_clim > vitot_min:
      vtot_target = vice_clim
    else:
      # Case when ithkn is turned off but also
      # this ignores hice = 0 in clim fields when ithkn is turned on
      vtot_target = np.max([vice_old, vice_new])

    ain_new, vin_new = adjust_thkncats_aice(ain_new, vin_new, vtot_target, \
                           hicat, dhi_min,  bnd_min=ain_min)
    
    ain_new = np.where(ain_new <= ain_min, 0., ain_new)
    vin_new = np.where(ain_new <= ain_min, 0., vin_new)

    assert np.sum(ain_new) < (1.+1e-12), f"Ice conc > 1 {np.sum(ain_new)} j={j0} i={i0}"

    # Check ice cats:
    hin_new = np.divide(vin_new, ain_new, out=np.zeros_like(vin_new), where=ain_new != 0)
    cat_missed, hcat_new = check_ithkn_cats(hicat, hin_new, ain_new)
    if cat_missed is not None:
      print(f"ipp={ipp} {aice_case} error in ice cats")
      raise Exception("Check ain_new, hin_new not in ice thkn cats")


    Tsfcn_new[:,j0,i0] = tsf_new
    vicen_new[:,j0,i0] = vin_new
    aicen_new[:,j0,i0] = ain_new

    # Update sice00?, qice00?
    for ilr in range(1, nilrs+1):
      varnum = f"{ilr:03d}"
      # Ice salinity by layers - compute S profile using BZ99 formulation:
      sice_lr  = sice_lr_cice4(ilr, nilrs, ain_new) 
      sice_old = sice[varnum][:,j0,i0]
      sice_new = np.where(sice_old < puny, sice_lr, sice_old)
      sice_new = np.where(ain_new < puny, 0., sice_new)

      sicen_new[varnum][:,j0,i0] = sice_new

      # ice enthalpy, use BL99
      # keep ice T below freezing T
      tice_max = sice_new*0 + Tfrz - 0.01
      qice_lr  = ice_enthalpy_BL99(tice_max, sice_new)
      qice_old = qice[varnum][:,j0,i0]
      qice_new = np.where(ain_old < puny, qice_lr, qice_old)
      qice_new = np.where(ain_new < puny, 0., qice_new)

      qicen_new[varnum][:,j0,i0] = qice_new

    # Check snow volume and ice_mean hsnow, add min snow if needed:
    vstot_old = np.sum(vsn_old)  # snow vol m3/m2_cell or cell mean snow thickn.
    if ai_new > puny:
      hsn_new = vstot_old / ai_new  # mean snow thickn over ice
      if hsn_new < hsnow_min:
        hsn_new = hsnow_min
        vstot_new = hsn_new * ai_new
      elif hsn_new > hsnow_max:
        hsn_new = hsnow_max
        vstot_new = hsn_new * ai_new
      else:
        vstot_new = vstot_old

      # Distribute evenly across cats:
      wts = ain_new/ai_new
      vsn_new = vstot_new * wts
    else:
      vsn_new = ain_new * 0.

    vsn_new = np.where(ain_new <= puny, 0., vsn_new)
    vsnon_new[:,j0,i0] = vsn_new

    # Update snow enthalpy: J/m3  
    # snow enthalpy should be: qsn_min <= qsn <= qsn_max
    # In theory, qsn_max = -rhos_Lfresh (latent heat of metling at 0C)
    # Make it a little lower to keep snow from melting right away
    #hsn_new = vsn_new / ain
    qsn = qsnon[:,j0,i0]        # enthalpy, J/kg < 0
    qsn_min = -rhos * Lfresh + (Tmin + 0.01) * cp_ice * rhos  # enth. of the coldest possible snow
    qsn_max = -rhos * Lfresh - 0.01 * cp_ice * rhos  # a little colder than 0C snow
    qsn_tsfc = -rhos * Lfresh + tsf_new * cp_ice * rhos # snow enth for surf. T
    qT0 = -Lfresh*rhos      # enth. of pure snow at 0C

    # Update new enthalpy of new snow:
    # Clip to min/max enthalpy, set to 0 where no snow:
    qsn_new = np.clip(qsn_tsfc, qsn_min, qsn_max)
    qsn_new = np.where(ain_new <= puny, 0., qsn_new)      # no ice
    qsn_new = np.where(vsn_new <= 0, 0., qsn_new)     # no snow
    qsnon_new[:,j0,i0] = qsn_new

  # Update data set:
  #ds_out = ds_out.assign(vsnon=vsnon_new, qsno001=qsnon_new)
  ds_out['vsnon'].values[:]   = vsnon_new
  ds_out['qsno001'].values[:] = qsnon_new
  ds_out['aicen'].values[:]   = aicen_new
  ds_out['vicen'].values[:]   = vicen_new
  ds_out['Tsfcn'].values[:]   = Tsfcn_new
  for ilr in range(1, nilrs+1):
    varnum = f"{ilr:03d}"
    ds_out[f"sice{varnum}"].values[:] = sicen_new[f"{varnum}"]
    ds_out[f"qice{varnum}"].values[:] = qicen_new[f"{varnum}"]

  # Sanity checking:
  assert ds_out["vsnon"].shape == vsnon_new.shape, "Check shape of vsnon "
  assert ds_out["qsno001"].shape == qsnon_new.shape, "Check shape of qsnon "

  check_print(aice, aicen_new, vicen_new, vsnon_new, ncat)

  return ds_out


def main():
  rest_date = None
  rest_hr = 0
  yrR = mmR = ddR = hrR = None
  yrN = mmN = ddN = hrN = None

  parser = argparse.ArgumentParser()
  parser.add_argument("--ithkn", type=int, default=1,
      help="insert ice thickn climatology, 0=no, 1=yes (default 1)", choices=[0,1])
  parser.add_argument("--rdate", help=f"restart date input file", required=True, type=int)
  parser.add_argument("--rhr", help=f"input file, restart hour = 0, ..., 23, default={rest_hr}", type=int)
  parser.add_argument("--rdate_out", help="output file, restart date", required=True, type=int)
  parser.add_argument("--rhr_out", help="output file, restart hour", required=True, type=int)
  parser.add_argument("--pth_in", help="input restart directory with original file", required=True, type=str)
  parser.add_argument("--flrst_in", help="rest file in", required=True, type=str)
  parser.add_argument("--pth_out", help="output restart directory where new file be dumped", required=True, type=str)
  parser.add_argument("--flrst_out", help="new rest file", required=True, type=str)
  parser.add_argument("--regn", help=f"where iconc incerted", 
                      choices=['north','south','global'], 
                      required=True, type=str)
  parser.add_argument("--fyaml",
    help=f"YAML with local directories, filenames, restart dates, default={default_yaml}, inside ./config",
    default=default_yaml,
    type=str)


  args = parser.parse_args()

  ins_thkn      = bool(args.ithkn)
  flrst_in      = args.flrst_in  
  flrst_out     = args.flrst_out 
  regn          = args.regn      
  rest_date     = args.rdate   
  rest_hr       = args.rhr       if args.rhr is not None else rest_hr
  rest_date_out = args.rdate_out 
  rest_hr_out   = args.rhr_out
  pth_in        = args.pth_in 
  pth_out       = args.pth_out
  fyaml         = CONFIG_DIR / args.fyaml

  print(f"Restart date input:  {rest_date}:{rest_hr}")
  print(f"Restart date output: {rest_date_out}:{rest_hr_out}")
  if ins_thkn:
    print("Insert NSDIC NRT ice concenatraion + ice thickness {fyaml} into CICE restart\n")
  else:
    print("Insert NSDIC NRT ice concenatraion, NO ice thickness\n")

  change_rest_time = (rest_date != rest_date_out) or (rest_hr != rest_hr_out)

  # Get date numbers:
  # Input restart file
  dnmbR = rdate2datenum(rest_date*100+rest_hr)  # restart day nmb
  if yrR is None:
    yrR,mmR,ddR,hrR = datevec(dnmbR, round_hrs=True)[:4]
  nsecR = hrR*3600

  # Dates of the output fields in the new restart:
  dnmbN = rdate2datenum(rest_date_out*100+rest_hr_out)
  if yrN is None:
    yrN, mmN, ddN, hrN = datevec(dnmbN, round_hrs=True)[:4]
  nsecN = hrN*3600
   
  if pth_in is None:
    raise RuntimeError("restart input path not provided")
  else:
    pthrest = pth_in

  # Output dir for new restart:
  if pth_out is None:
    raise RuntimeError("restart output path not provided")
  else:
    pthrest_out = pth_out

  with open(fyaml) as ff:
    config_rest = yaml.safe_load(ff)

  def expand_env_vars(d):
      for k, v in d.items():
          if isinstance(v, dict):
              expand_env_vars(v)
          elif isinstance(v, str):
              d[k] = os.path.expandvars(v)

  expand_env_vars(config_rest)

  # Get MOM6 grid
  pthgrid    = config_rest["grid_topo"]["pthgrid"]
  dfgrid_mom = os.path.join(pthgrid, "ocean_hgrid.1440x1080.nc")
  dftopo_mom = os.path.join(pthgrid, "ocean_topog.1440x1080.nc")

  with xarray.open_dataset(dftopo_mom) as dstopo:
    HH = dstopo['depth'].data.squeeze()

  HH = np.where(HH < 1.e-20, np.nan, HH)
  HH = -HH
  HH = np.where(np.isnan(HH), 1., HH)
  jdm, idm = HH.shape

  LON, LAT = read_mom6grid(dfgrid_mom, grdpnt='hgrid')

  print(f"old restart: {yrR}/{mmR:02d}/{ddR:02d}:{hrR:02d}")
  print(f"new restart: {yrN}/{mmN:02d}/{ddN:02d}:{hrN:02d}")

  ds_out = None
  if regn == 'global':
    for regn_tmp in (['north','south']):
      print(f"  Processing region: {regn_tmp}")
      ds_out = insert_iconc_ithkn(regn_tmp, ds_out, config_rest, HH, LAT, LON, dnmbR, dnmbN,\
                       ins_thkn, pthrest, flrst_in)
  else:
    print(f"  Processing region: {regn}")
    ds_out = insert_iconc_ithkn(regn, ds_out, config_rest, HH, LAT, LON, dnmbR, dnmbN,\
                       ins_thkn, pthrest, flrst_in)

  # Attributes:
  if ins_thkn:
    title_str = f"CICE6 restart with inserted ice conc NSIDC NRT {rest_date_out} and ice thickness from {fyaml}"
  else:
    title_str = f"CICE6 restart with inserted ice concentration from NSIDC NRT {rest_date_out} "

  from datetime import datetime
  istep1_val = ds_out.attrs.get('istep1', None)
  ds_out.attrs.update({
      "title": title_str,
      "source": "insert_iconc_ithkn_cice6rest_global.py",
      "istep1": np.int32(istep1_val) if istep1_val is not None else np.int32(0), 
      "myear": np.int32(yrN),
      "mmonth": np.int32(mmN),
      "mday": np.int32(ddN),
      "msec": np.int32(nsecN),
      "history1": f"Modified {flrst_in} {datetime.now().isoformat()}",
      "region": regn,
  })

  # Save:
  if flrst_out is None:
    if ins_thkn:
      flrst_out = f"cice_model.res.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.iconc_thkn.nc"
    else:
      flrst_out = f"cice_model.res.{yrN}{mmN:02d}{ddN:02d}.{hrN:02d}.iconc.nc"

  dflrst_out = os.path.join(pthrest_out, flrst_out)
  print(f"Saving CICE restart --> {dflrst_out}")
  ds_out.to_netcdf(dflrst_out, encoding={var: {'_FillValue': None} for var in ds_out.data_vars}, format='NETCDF3_64BIT')
  ds_out.close()


if __name__ == "__main__":
  main()


