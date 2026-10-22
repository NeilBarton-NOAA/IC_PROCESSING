"""
  MOM6 utilities
  reading grid
  Dmitry Dukhovskoy, NOAA OAR PSL 
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import xarray
from netCDF4 import Dataset as ncFile

def read_mom6grid(fgrid, grdpnt='hgrid', grid='nonsymmetr'):
  """
    Read MOM6 grid
    coordinates on "supergrid" (half grid points of actual grid)
    default - coordinates of h-points (scalar) are derived
    options: 'qgrid' - coordinated in the upper-right corner
                      q-points, vorticity points

   https://mom6.readthedocs.io/en/main/api/generated/pages/Discrete_Grids.html
 
   q(i-1,j)   ------------  v(i,j) --------------  q(i,j)
     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  

   u(i-1,j)               * h(i,j)                u(i,j)

     |                                              |  
     |                                              |  
     |                                              |  
     |                                              |  
   q(i-1,j-1) ------------  v(i,j-1) ------------  q(i,j-1)


  grid nonsymmetr: All arrays are declared with the same shape (isd:ied,jsd:jed)
                 There are fewer staggered variables to the south-west of the 
                 computational domain. An operator applied at h-point locations 
                 involving u- or v- point data can not have as wide a stencil on 
                 the south-west side of the processor domain as it can on 
                 the north-east side.

  grid symmetr:   Arrays have different shapes depending on their staggering 
                  location on the Arakawa C grid. 

  """

  print(f"Reading MOM6 grid for: {grid} grdpnt={grdpnt}")
  print('Grid: ' + fgrid)
  nc  = ncFile(fgrid,'r')
# Read MOM supergrid:
  XX  = nc.variables['x'][:].data
  YY  = nc.variables['y'][:].data
  mm  = XX.shape[0]
  nn  = XX.shape[1]


  if grid == 'nonsymmetr':
    if grdpnt == 'qgrid':
      LON = XX[2::2, 2::2]
      LAT = YY[2::2, 2::2]
    elif grdpnt == 'hgrid':
      LON = XX[1::2, 1::2]
      LAT = YY[1::2, 1::2]
    elif grdpnt == 'ugrid':
      LON = XX[1::2, 2::2]
      LAT = YY[1::2, 2::2]
    elif grdpnt == 'vgrid':
      LON = XX[2::2, 1::2]
      LAT = YY[2::2, 1::2]
  elif grid == 'symmetr':
    if grdpnt == 'qgrid':
      LON = XX[0::2, 0::2]
      LAT = YY[0::2, 0::2]
    elif grdpnt == 'hgrid':
      LON = XX[1::2, 1::2]
      LAT = YY[1::2, 1::2]
    elif grdpnt == 'ugrid':
      LON = XX[1::2, 0::2]
      LAT = YY[1::2, 0::2]
    elif grdpnt == 'vgrid':
      LON = XX[0::2, 1::2]
      LAT = YY[0::2, 1::2]
  elif grid == 'supergrid':
# supergrid is the same for all pnts, 
    LON = XX
    LAT = YY

  jdm = LAT.shape[0]
  idm = LAT.shape[1]

  print(f" MOM6 {grdpnt} {grid} grid dim (J,I): {jdm} x {idm}")
  LON = np.where(LON < -180., LON+360., LON)

  return LON, LAT

def read_mom6angle(fgrid, grdpnt='hgrid', grid='nonsymmetr'):
  """
    Read angle the grid makes with the true East direction (lat line)
    radians
  """
  print(f"Reading MOM6 grid angle (radians) for: {grid} grdpnt={grdpnt}")
  print('Grid: ' + fgrid)
  nc  = ncFile(fgrid,'r')
# Read MOM supergrid:
  angle = nc.variables['angle_dx'][:].data
  mm  = angle.shape[0]
  nn  = angle.shape[1]

  if grid == 'nonsymmetr':
    if grdpnt == 'qgrid':
      alpha = angle[2::2, 2::2]
    elif grdpnt == 'hgrid':
      alpha = angle[1::2, 1::2]
    elif grdpnt == 'ugrid':
      alpha = angle[1::2, 2::2]
    elif grdpnt == 'vgrid':
      alpha = angle[2::2, 1::2]
  elif grid == 'symmetr':
    if grdpnt == 'qgrid':
      alpha = angle[0::2, 0::2]
    elif grdpnt == 'hgrid':
      alpha = angle[1::2, 1::2]
    elif grdpnt == 'ugrid':
      alpha = angle[1::2, 0::2]
    elif grdpnt == 'vgrid':
      alpha = angle[0::2, 1::2]

  jdm = alpha.shape[0]
  idm = alpha.shape[1]

  print(f" MOM6 dx angle {grdpnt} {grid} grid dim (J,I): {jdm} x {idm}")

  return  alpha

def read_mom6depth(ftopo, varnm='depth', f_negate=True):
  """
    Read MOM6 depths at h-pnts
  """
  print(f'Reading MOM6 depths {ftopo}')
  with xarray.open_dataset(ftopo) as dnc:
    HH = dnc[varnm].data.squeeze()

# Convert depth to negatives:
# Land > 0
  if f_negate:
    HH  = np.where(HH > 0, -HH, 100.)

  return HH

def read_mom6lmask(ftopo):
  """
    Read MOM6 land mask
    1 = ocean, 0 = land
  """
  print('Reading MOM6 land mask ' + ftopo)
  nc   = ncFile(ftopo,'r')
  Lmsk = nc.variables['wet'][:].data

  return Lmsk


def read_mom6(foutp, fld, rLayer=-1, fnan=True, finfo=True):
  """
    Read 2D or 3D field from MOM6 output
    Specify Lyaer #  to read for 3D fields (from 1, ..., kdm)
    otherwise all layers are read in 3D array
    Replace out of range data with nan's - default

  """
  print('Reading MOM6 {0}: {1}'.format(fld, foutp))

  huge = 1.e18
  nc   = ncFile(foutp, 'r')

# Check variable dimension
# Assumed: AA(Time, dim1, dim2, optional: dim3)
  ndim = nc.variables[fld].ndim - 1 # Time excluded
  
# Check that time dim=1:
  tdim = nc.variables['Time'][:].data.shape[0]
  if tdim != 1:
    raise Exception('in netcdf output  Time dimensions is not 1: {0}'.\
                    format(tdim))

  if ndim == 2:
    FF = nc.variables[fld][:].squeeze().data
  elif ndim == 3:
    if rLayer <= 0:
      FF  = nc.variables[fld][:].squeeze().data
    else:
      klr = rLayer-1
      FF  = nc.variables[fld][:,klr,:,:].squeeze().data
  else:
    print('Output field dim should be 2 or 3')
    raise Exception('Check {0} dimension = {1}'.format(fld,ndim))

  AA = np.where(FF >= huge, np.nan, FF)
  if finfo:
    print('{0} lr {3} min/max: {1} / {2}'.\
         format(fld,np.nanmin(AA),np.nanmax(AA),rLayer))

  if fnan:
    FF = AA

  return FF

def mom_dim(flin):
  """
    Get horiz/vert dimensions for mom6
  """  
  nc   = ncFile(flin, 'r')
  xh = nc.variables['xh'][:].data.shape[0]
  yh = nc.variables['yh'][:].data.shape[0]
  zl = nc.variables['zl'][:].data.shape[0]

  return xh, yh, zl

def zz_zm_fromDP(dH, ssh, f_intrp=False, f_btm=True, f_ssh = True, \
                 finfo=True, eps0=1.e-5):
  """
    Calculate ZZ, ZM from layer thkcness (dH) 
    ZZ - interface depths,
    ZM - mid cell depths
    Note that in MOM sum(dH) = water column height that includes SSH
    Therefore ZZ[0] = ssh

    f_btm = true - nan below bottom & land
            false - no nans, all ZZ=zbtm or 0
    f_intrp = True: for vertical interpolation, land/bottom make
              not nans and not equal depths, i.e. depth continue
              increasing at eps0 rate below bottom to keep
              monotonicity of zz, zm
    f_ssh = False - assumes ZZ[0] = 0, i.e. ssh effect should be removed
                     from lr. thicknesses, see:
                     remove_ssh_lrthk
    f_interp overrides f_btm flag: makes it False
  """
  if f_intrp:
    f_btm = False

  print('Deriving ZZ, ZM from dH ...')
  ll = dH.shape[0]
  mm = dH.shape[1]
  nn = dH.shape[2]

  ZZ = np.zeros((ll+1,mm,nn))
  ZM = np.zeros((ll,mm,nn))
  ssh = np.where(np.isnan(ssh),0.,ssh)
  
  if f_ssh:
    ZZ[0,:,:] = ssh
  else:
    ZZ[0,:,:] = 0.

  if f_btm:
    dH = np.where(np.isnan(dH),0., dH)
  else:
    dH = np.where(np.isnan(dH),eps0, dH)

  if f_intrp:
    dH = np.where(dH < eps0, eps0, dH)

  for kk in range(ll):
    ZZ[kk+1,:,:] = ZZ[kk,:,:]-dH[kk,:,:]

    if f_btm:
      [JJ,II] = np.where(dH[kk,:,:] <= eps0)
      ZZ[kk+1,JJ,II] = np.nan
      if kk==0:
        ZZ[kk,JJ,II] = np.nan
    elif not f_btm and not f_intrp:
      [JJ,II] = np.where(dH[kk,:,:] <= eps0)
      ZZ[kk+1,JJ,II] = ZZ[kk,JJ,II]

# Depths of the middle of the grid cells:
    ZM[kk,:,:] = 0.5*(ZZ[kk+1,:,:]+ZZ[kk,:,:])
    if finfo:
      print(' kk={0} min/max ZZ = {1:5.1f}/{2:5.1f}'.\
           format(kk+1,np.nanmin(ZZ[kk,:,:]),np.nanmax(ZZ[kk,:,:])))

  return ZZ, ZM

def zm2zz (ZM):
  """
    Derive layer interface depths from layer mid-points
    1st dimension is depth for ZM
  """
  ndim = len(ZM.shape)
  if ndim > 3:
    raise Exception(f"Input dZ dim={ndim}, cannot be > 3")

  dim2 = 1
  dim3 = 1
  if ndim == 1:
    dim1 = ZM.shape[0]
  elif ndim == 2:
    dim1, dim2 = ZM.shape
  elif ndim == 3:
    dim1, dim2, dim3 = ZM.shape

  ZZ = np.zeros((dim1+1,dim2,dim3)).squeeze()
  zz_negate = np.min(ZM) < 0.  

  for kk in range(dim1):
    if ndim == 1:
      dz = abs(ZM[kk] - ZZ[kk])
      if zz_negate: dz = -dz
      ZZ[kk+1] = ZM[kk] + dz
    else:
      dz = abs(ZM[kk,:] - ZZ[kk,:])
      if zz_negate: dz = -dz
      ZZ[kk+1,:] = ZM[kk,:] + dz
    
  return ZZ  

def zz2zm (ZZ):
  """
    Derive layer mid- depths from layer interface depths
    1st dimension is depth for ZZ
  """
  ndim = len(ZZ.shape)
  if ndim > 3:
    raise Exception(f"Input dZ dim={ndim}, cannot be > 3")

  dim2 = 1
  dim3 = 1
  if ndim == 1:
    dim1 = ZZ.shape[0]
  elif ndim == 2:
    dim1, dim2 = ZZ.shape
  elif ndim == 3:
    dim1, dim2, dim3 = ZZ.shape

  ZM = np.zeros((dim1-1,dim2,dim3)).squeeze()
  zm_negate = np.min(ZZ) < 0.

  for kk in range(dim1-1):
    if ndim == 1:
      dz = abs(ZZ[kk+1] - ZZ[kk])
      if zm_negate: dz = -dz
#      print(f'kk={kk} dz={dz:5.2f} ZZ={ZZ[kk]}')
      ZM[kk] = ZZ[kk] + 0.5*dz
    else:
      dz = abs(ZZ[kk+1,:] - ZZ[kk,:])
      if zm_negate: dz = -dz
      ZM[kk,:] = ZZ[kk,:] + 0.5*dz

  return ZM

def zz_zm_fromDZ(dZ, depth_negate = True):
  """
    Simple algorithm for getting layer interface depths (zz) and mid-point depths (zm)
    from layer thicknesses (dZ)
    in dZ 1st dimension is depth 
    No bottom / ssh 
  """
  ndim = len(dZ.shape)
  if ndim > 3:
    raise Exception(f"Input dZ dim={ndim}, cannot be > 3")

  dim2 = 1
  dim3 = 1
  if ndim == 1:
    dim1 = dZ.shape[0]
  elif ndim == 2:
    dim1, dim2 = dZ.shape
  elif ndim == 3:
    dim1, dim2, dim3 = dZ.shape

  ZZ = np.zeros((dim1+1,dim2,dim3)).squeeze()
  ZM = np.zeros((dim1,dim2,dim3)).squeeze()
  if depth_negate: dZ = -(abs(dZ))

  if ndim == 1:
    for kk in range(dim1): 
      ZZ[kk+1] = ZZ[kk] + dZ[kk]
      ZM[kk]   = 0.5*(ZZ[kk+1] + ZZ[kk])
  else:
    for kk in range(dim1): 
      ZZ[kk+1,:] = ZZ[kk,:] + dZ[kk,:]
      ZM[kk,:]   = 0.5*(ZZ[kk+1,:] + ZZ[kk,:])

  return ZZ, ZM
    

def get_zz_zm(fina, f_btm=True, **kwargs):
  """
    Derive layer depths: ZZ - interface depths,
    ZM - mid cell depths
    Note that in MOM sum(dH) = water column height that includes SSH
    Therefore ZZ[0] = ssh
    f_btm = true - ZZ=nan below bottom

    can specify the names for layer thickness, ssh variables
    if these are different from 'h' and 'SSH' as kwargs, e.g."
    [sshvar="ssh", lrthkvar="hlr"]
  """
  print(' Deriving ZZ, ZM from '+fina)
  huge = 1.e18
  nc   = ncFile(fina, 'r')

  sshvar   = 'SSH'
  lrthkvar = 'h'
  for key, value in kwargs.items():
    if key == 'sshvar':
      sshvar = value
    elif key == 'lrthkvar':
      lrthkvar = value 

# Check variable dimension
# Assumed: AA(Time, dim1, dim2, optional: dim3)
  ndim = nc.variables['h'].ndim - 1 # Time excluded

# Check that time dim=1:
  try:
    tdim = nc.variables['Time'][:].data.shape[0]
  except:
    try:
      tdim = nc.variables['time'][:].data.shape[0]
    except:
      print(f"variable Time/time not found in {fina}")

  if tdim != 1:
    raise Exception('in netcdf output  Time dimensions is not 1: {0}'.\
                    format(tdim))

# Read layer thicknesses:
  dH     = nc.variables[lrthkvar][:].squeeze().data
  dH     = np.where(dH >= huge, np.nan, dH)
  ssh    = nc.variables[sshvar][:].squeeze().data
  ssh    = np.where(ssh >= huge, np.nan, ssh)
  ZZ, ZM = zz_zm_fromDP(dH, ssh, f_btm=f_btm)  

  return ZZ, ZM

def remove_ssh_lrthk(ssh, HH, dH):
  """
  Remove ssh from layer thicknesses such that
  sum_over_z(dH[:,j,j])-HH[j,i] = 0 otherwise sum(dH[:,j,i])-HH[j,i] = ssh

  dH is 3D array of layer thicknesses
  HH - botom bathymetry 2D
  ssh - sea surface height 2D

  SSH correction to thickness:
  See HYCOM-tools: archv2mom6res.f
  qq  = (ssh0 + abs(hbH))/abs(hbH)
  dh0 = dh0*qq 

  """
  qqm  = (ssh + abs(HH))/abs(HH)
  kdim = dH.shape[0]
  dH0  = dH.copy()
  for ik in range(kdim):
    aa = dH[ik,:,:]
    aa = aa/qqm     # remove ssh correction from layers
    dH0[ik,:,:] = aa

  return dH0

def create_time_array(date1, date2, dday, date_mat=False):
  """
    Create time array for plotting fields
    if date_mat: date1 and date2 are in mat format
    otherwise: dates = [year, month, day]

  """
  import mod_time as mtime
  importlib.reload(mtime)

  if isinstance(date1, float) or isinstance(date1, int):
    date_mat = True

  if not date_mat:
#    yr1, mo1, mday1 = date1 
#    yr2, mo2, mday2 = date2  
#
    dnmb1 = mtime.datenum(date1)
    dnmb2 = mtime.datenum(date2)
  else:
    dnmb1 = date1
    dnmb2 = date2

  ldays = np.arange(dnmb1,dnmb2+dday/100.,dday)   
  nrec  = ldays.shape[0]
  TPLT  = np.zeros(nrec, dtype=[('dnmb', float),
                                ('date', int, (4,)),
                                ('yrday', float)])

  for irc in range(nrec):
    dnmb = ldays[irc]
    DV   = mtime.datevec(dnmb)
    _, jday = mtime.dnmb2jday(dnmb)

    TPLT['dnmb'][irc]   = dnmb
    TPLT['yrday'][irc]  = jday
    TPLT['date'][irc,:] = DV[0:4]

  print('Creating TPLT Time array, start: ' + \
        '{0}/{1}/{2} {3}hr, '.format(TPLT['date'][0,0], \
        TPLT['date'][0,1], TPLT['date'][0,2], TPLT['date'][0,3]) + \
        'end: {0}/{1}/{2} {3}hr'.\
      format(TPLT['date'][-1,0], TPLT['date'][-1,1], TPLT['date'][-1,2],\
      TPLT['date'][-1,3]))

  return TPLT  


def dx_dy(LON,LAT):
  """
    Find horizontal grid spacing from LON,LAT 2D arrays 
    hycom grid
  """
  import mod_misc1 as mmsc1

  IDM = LON.shape[1]
  JDM = LON.shape[0]
  print('Calculating DX, DY for HYCOM idm={0}, jdm={1}'.format(IDM,JDM))

  DX = np.zeros((JDM,IDM))
  DY = np.zeros((JDM,IDM))

  for ii in range(IDM-1):
    LT1 = LAT[:,ii]
    LT2 = LAT[:,ii+1]
    LN1 = LON[:,ii]
    LN2 = LON[:,ii+1]
    dx  = mmsc1.dist_sphcrd(LT1,LN1,LT2,LN2)
    DX[:,ii] = dx

  DX[:,IDM-1] = dx

  for jj in range(JDM-1):
    LT1 = LAT[jj,:]
    LT2 = LAT[jj+1,:]
    LN1 = LON[jj,:]
    LN2 = LON[jj+1,:]
    dy  = mmsc1.dist_sphcrd(LT1,LN1,LT2,LN2)
    DY[jj,:] = dy

  DY[JDM-1,:] = dy

  print('Min/max DX, m = {0:12.4f} / {1:12.4f}'.format(np.min(DX),np.max(DX)))
  print('Min/max DY, m = {0:12.4f} / {1:12.4f}'.format(np.min(DY),np.max(DY)))

  return DX, DY

 
def vol_transp_2Dsection(LSgm, ZZ, UV):
  """
    Compute Vol transport along a straight line (EW or SN orientation)
    Note this algorithm works for zigzaging section
    as long as ZZ, UV are collocated in the center of the segments
    of the contour line
   
    Use midpoint quadrature
    All inputs are 2D fields (ZZ can be 1D):
    LSgm - segment lengths (grid cell dx)
    ZZ   - interface depths, can be 1D array
    UV   - normal U component 
  """
  klv, nsgm  = UV.shape
  if len(np.shape(ZZ)) == 1:
    Z1d = ZZ.copy()
    for ik in range(nsgm-1):
      ZZ = np.column_stack((ZZ, Z1d))
  dZ   = abs(np.diff(ZZ, axis=0))

  Ctrp = np.zeros((klv,nsgm))
  for kk in range(klv):
    Ctrp[kk,:] = UV[kk,:]*dZ[kk,:]*LSgm

  volTr_1D = np.nansum(Ctrp, axis=0)

  return volTr_1D

def fill_bottom(A2d, ZZ, Hbtm, fill_land=True):
  """
    Fill bottom values for plotting
    ZZ - layer interface depths 1D or 2D
    fill_land - fill land or not with closest values
  """
  A2df = A2d.copy()
  
  if len(ZZ.shape) == 2:
    f2d = True
  else:
    f2d = False
    z0 = ZZ

  kdm = A2d.shape[0]
  idm = A2d.shape[1]
  iOc = np.where(Hbtm < 0.)[0]

  for ii in range(idm):
    hb0 = Hbtm[ii]
    if hb0 >= 0.:
      if fill_land:
        dii        = abs(iOc - ii)
        ix         = np.where(dii == np.min(dii))[0][0]
        iocn       = iOc[ix]
        A2df[:,ii] = A2d[:,iocn] 
        continue
      else:
        continue 

    if f2d:
      z0 = ZZ[:,ii]

    dZb = ZZ-hb0
#    print(f'ii={ii} minZZ={np.min(ZZ):4.1f} hb={hb0:4.1f}')
# If bottom deeper than the deepest interface:
    if abs(ZZ[-1]) < abs(hb0): continue
    izb = min(np.where(dZb <= 1.e-3)[0])

    A2df[izb:,ii] = A2d[izb-1,ii]
    
  return A2df

def bottom2nan(A2d, ZZ, Hbtm):
  """
    Fill bottom/land values with nans for plotting
    ZZ - layer interface depths 1D or 2D
  """
  A2df = A2d.copy()
  
  if len(ZZ.shape) == 2:
    f2d = True
  else:
    f2d = False
    z0 = ZZ

  kdm = A2d.shape[0]
  idm = A2d.shape[1]

  for ii in range(idm):
    hb0 = Hbtm[ii]
    if hb0 >= 0.:
     A2df[:,ii] = np.nan

    if f2d:
      z0 = ZZ[:,ii]

    dZb = ZZ-hb0
    izb = min(np.where(dZb <= 1.e-3)[0])
    A2df[izb:,ii] = np.nan
    
  return A2df

def collocateU2H(A2d, grid_shape, f_land0 = True):
  """
    Collocate variables from MOM u-point to H-point
    input 2D field (1 layer) at p-point
    Collocation is done "as is", i.e.
    land values will be brought into H-point as is
    If these are nans or some filled values - double check
    f_land0 = Replace land values (nans) with 0 - ok for u/v but 
              not for T/S

    grid_shape = symmetr/ nonsymmetr

    MOM6 grid:

     Q(i-1,j) V(i,j)     Q(i,j)
     * -------|----------*
     |                   |
     |        P,T,S,H    |
     |        * i,j      - U(i,j)
     |                   |
     |                   |
     * ------------------*
   Q(i-1,j-1)        

  """
  inan = np.where(np.isnan(A2d))[0]
  if len(inan) > 0 and not f_land0:
    print('collocateU2H: WARNING - nans in U, may impact near-coast U values, use f_land0')

  lsymmetr = False
  if grid_shape[0:7] == 'symmetr':
    lsymmetr = True

  mm  = A2d.shape[0]
  nn  = A2d.shape[1]
  mhp = mm
  nhp = nn 
  if lsymmetr: nhp=nn-1

  if f_land0:
    A2d  = np.where(np.isnan(A2d), 0., A2d)
    A2d  = np.where(A2d > 1.e10, 0., A2d)
  A2dP = np.zeros((mhp, nhp))

  if lsymmetr:
    for ii in range(1,nn):
      u1 = A2d[:,ii]
      u2 = A2d[:,ii-1]
      uP         = 0.5*(u1 + u2)
      A2dP[:,ii-1] = uP
  else:
    for ii in range(nn):
      u1 = A2d[:,ii]
      if ii > 0:
        u2 = A2d[:,ii-1]
      else:
        u2 = A2d[:,ii]    # no periodic open boundaries assumed
      uP         = 0.5*(u1 + u2)
      A2dP[:,ii] = uP

  return A2dP

def collocateV2H(A2d, grid_shape, f_land0 = True):
  """
    Collocate variables from MOM v-point to H-point
    input 2D field (1 layer) at p-point
    Collocation is done "as is", i.e.
    land values will be brought into H-point as is
    If these are nans or some filled values - double checl
    f_land0 = Replace land values (nans) with 0 - good for u/v but 
              not for T/S

    grid_shape = symmetr/ nonsymmetr

    HYCOM grid:

     Q(i-1,j) V(i,j)     Q(i,j)
     * -------|----------*
     |                   |
     |        P,T,S,H    |
     |        * i,j      - U(i,j)
     |                   |
     |                   |
     * ------------------*
   Q(i-1,j-1)        

  """
  inan = np.where(np.isnan(A2d))[0]
  if len(inan) > 0 and not f_land0:
    print('collocateV2H: WARNING - nans in V, may impact near-coast V values, use f_land0')

  lsymmetr = False
  if grid_shape[0:7] == 'symmetr':
    lsymmetr = True
    
  mm  = A2d.shape[0]
  nn  = A2d.shape[1]
  mhp = mm
  nhp = nn
  if lsymmetr: mhp=mm-1

  if f_land0:
    A2d  = np.where(np.isnan(A2d), 0., A2d)
    A2d  = np.where(A2d > 1.e10, 0., A2d)
  A2dP = np.zeros((mhp, nhp))

  if lsymmetr:
    for jj in range(1,mm):
      v1 = A2d[jj,:]
      v2 = A2d[jj-1,:]
      vP         = 0.5*(v1 + v2)
      A2dP[jj-1,:] = vP
  else:
    for jj in range(mm):
      v1 = A2d[jj,:]
      if jj > 0:
        v2 = A2d[jj-1,:]
      else:
        v2 = A2d[jj,:]    # no periodic open bovndaries assvmed
      vP         = 0.5*(v1 + v2)
      A2dP[jj,:] = vP

  return A2dP

def fill_land3d(A3d, vert2d=False,  **kwargs):
  """
    Fill NaNs with closest values
    2D or 3D arrays allowed
    by default 2D is horizontal field and no bottom filled values performed
    for vertical 2D sections make vert2d True

    First, the 1st layer is filled - land values interpolated from the closest ocean points
    The values are smoothed to avoid jumps near the ocean land points
    Next (for 3D arrays) - below 1st layer: values are filled with the above lr value
 
    For quick fill, specify quick_fill=Value, all Nans will be filled with Value

    Note: bottom values at the ocean grid points are not checked
          but nans will be filled for quick_fill with a specified value
          for not quick_fill - with values above nans

    Filled values can be smoothed - boxfltr = squared box size (should be = n^2, n = odd 3, 5, 7,...)

    land_mask: value:  Replace land points masked as some huge=land_mask value with nans

    verbosity = 1: print info about dim of the array being filled, =0 - silent

    Usage: Afilled = mod_mom6.fill_land3d(A3d, [quick_fill=1.e22, boxfltr=25]) 
  """
  import mod_interp1D as minterp
  import mod_misc1 as mmisc
  from scipy import interpolate

  qfill = False
  bxflt = False
#  bxflt = True
#  nbx = 
  huge = 1.e20   # possible land values if not nans
  huge2nan = False
  verbosity = 1
  sinfo = ''
  for key, value in kwargs.items():
    if key == 'quick_fill': 
      vfill = value
      qfill = True
    if key == 'boxfltr':
      nbx = value
      bxflt = True
    if key == 'land_mask':
      huge = value
      huge2nan = True
    if key == 'verb':
      verbosity = value
    if key == 'sinfo':
      sinfo = value

  # Check if there are any land points, should be nan's
  aa = np.max(abs(A3d))
  if ~np.isnan(aa) and aa < huge:
    print(f'No land values found, no land filling required max val={aa}')
    return A3d
  elif aa >= huge:
    if huge2nan:
      A3d = np.where(A3d >= huge, np.nan, A3d)
    else:
      print(f'WARNING: land values should be nans, found max val={aa}, is this land?')
      print(f'WARNING: use keyarg land_mask={aa} to change to nan, no land filling performed ...')
      return A3d

  adim = A3d.shape
  ndim = len(adim)
  if ndim > 3 or ndim < 2:
    raise Exception('Array should be 2 or 3D')
  if verbosity > 0:
    print(f'{sinfo} Filling land, {ndim}D array')

  f2d = False
  if ndim == 3:
    kdim = adim[0]
    jdim = adim[1]
    idim = adim[2]
  else:
    f2d = True
    kdim = 0
    jdim = adim[0]
    idim = adim[1]
 
  if qfill:
    print(f'Quick land fill option with val={vfill}')
    A3d = np.where(np.isnan(A3d), vfill, A3d)
    return A3d

# First, fill land in surface layer: 
  if f2d:
    A2d = A3d
  else:
    A2d = A3d[0,:,:].squeeze()

  J,I = np.where(~np.isnan(A2d))
  ny, nx = A2d.shape
  X = np.arange(0,nx)
  Y = np.arange(0,ny)
  XX, YY = np.meshgrid(X,Y)
  A2df = interpolate.griddata((I,J), A2d[J,I], (XX,YY), method='linear')

# interpolation leaves nans in the corners
  Nnan = len(np.where(np.isnan(A2df))[0])
  if Nnan > 0:
    for jj in range(jdim):
      imss = np.where(np.isnan(A2df[jj,:]))[0]
      nmss = len(imss)
      if nmss == 0: continue
      a1d  = A2df[jj,:] 
      i1d  = np.where(~np.isnan(a1d))[0]
      v1d  = A2df[jj,i1d]
  # Add endpoints for interpolation:
  #    print(f'jj={jj} len i1d={len(i1d)}')
      if len(i1d) == 0:
  # Case all nans - wall OB or 2D vertical section below bottom
  # try above row or overall mean
        v1d = np.zeros((2))
        i1d = np.zeros((2))
        i1d[0] = 0
        i1d[1] = idim-1
        if jj > 0:
          v1d[0] = A2df[jj-1,0]
          v1d[1] = A2df[jj-1,-1]
        else:
          v1d[0] = np.nanmean(A2df)
          v1d[1] = np.nanmean(A2df)
   
      if not i1d[0] == 0:
        v1d = np.insert(v1d,0,v1d[0])
        i1d = np.insert(i1d,0,0)

      if i1d[-1] < idim-1:
        v1d = np.append(v1d, v1d[-1])
        i1d = np.append(i1d, idim)

      for I0 in range(nmss):
        ixx = imss[I0]
        xF  = minterp.pcws_lagr1(i1d,v1d,ixx)
        a1d[ixx] = xF

  # Check that ocean points have not been impacted:
      dltA = np.nanmax(np.abs(A2df[jj,:]-a1d))
      if dltA > 1.e-9:
        raise Excpetion(f"Ocean point contaminated during land filling j={jj} dltA={dltA}")

      A2df[jj,:] = a1d

  if bxflt:
    JN,IN = np.where(~np.isnan(A2d))
    A2df = mmisc.box_fltr(A2df, nbx=nbx)
    A2df[JN,IN] = A2d[JN,IN]

  if f2d: 
    return A2df

  A3df = A3d.copy()
  A3df[0,:,:] = A2df

# Fill deep layers:
#  print('Start filling deep laeyrs')
  for kk in range(1,kdim):
    A2d = A3df[kk,:,:].squeeze()
    JN,IN = np.where(np.isnan(A2d))
    nall = len(JN)
#    print(f'Layer {kk+1} found nans={nall}')
    if nall == 0: continue
    A3df[kk,JN,IN] = A3df[kk-1,JN,IN]

  return A3df


def deallocateP2UV(A2d, grid_shape, f_component, f_ignore=True, **kwargs):
  """
    Deallocate U or V component from P point putting it back to V-grid
    Best we can do - use interpolation hoping to get close enough
    to the true value at v-point
    Specify via keyword arguments polynomial: 2nd or 3rd 
    Default - cubic polynomial
    f_component = 'v-vel' or 'u-vel'

   
    x1 v-value at v-pnt (i1) 
    *


         o p-point      * x3 (i1+2)
        value xp1 
        (i1+1/2)   o p-point (i1+1.5)
                      
              * x2 (i1+1) 

    Given values at p-points - need to estimate values (x1, x2, ..) at v-pnts 
    A2d - 2D array with v at p-points
 
    Use polynom. interpolation 
     
   ==== NOT FINISHED =====


  """
  print('UNFINISHED ....')
  return

  f_interp = 'cubic'
  for key, value in kwargs.items():
    if key == "interp": f_interp = value

  if f_interp == 'cubic':
    nnd = 4
  else:
    nnd = 3

  djm = int(np.floor((nnd-1)/2))
  djp = nnd-djm-1

# Interpolate  along the rows
  if f_component == 'u-vel':
    A2d = A2d.transpose()

  idim = A2d.shape[1]
  jdim = A2d.shape[0]
  A2di = np.zeros((jdim, idim))

  for jj in range(jdim):
    if jj == jdim-1:
      A2di[jj,-1] = A2d[jj,-1]
      
    J = np.arange(0,nnd)+jj-djm
    if min(J) < 0: J = J-min(J)  
# at the "right" boundary use piecewise Lagr cardinal basis
# on an interval left of xi x=[x(i-nnd), xi], order of nodes
# is inverse, so that x(i-1) = 
    if max(J) >= jdm:
      djj = max(J)-jdm+1
      J = np.flip(J-djj)
  
#     P2d = minterp.lagr_polynom2D
  if f_component == 'u-vel': A2di = A2di.transpose()

  return


