"""
  Utilities 
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import pdb
import importlib
import struct
import pickle
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
#from mpl_toolkits.basemap import Basemap, shiftgrid

sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/hycom_utils')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython/draw_map')
sys.path.append('/home/Dmitry.Dukhovskoy/python/MyPython')
import mod_read_hycom
#importlib.reload(mod_read_hycom)
from mod_read_hycom import read_grid_topo, read_hycom, read_topo
from mod_read_hycom import zz_zm_fromDP
from mod_utils_fig import bottom_text

def ncoda_depths(zneg=True):
  """
    NCODA interface and mid-point depths
    zi = from hycom/ncoda_archv_inc/zi.txt
  """
  ZI = np.array([      
      0.00,
      1.00,
      3.00,
      5.00,
     11.00,
     21.00,
     35.00,
     53.00,
     67.00,
     85.00,
     99.00,
    117.00,
    131.00,
    149.00,
    171.00,
    189.00,
    211.00,
    229.00,
    251.00,
    269.00,
    291.00,
    309.00,
    371.00,
    389.00,
    451.00,
    469.00,
    531.00,
    589.00,
    651.00,
    709.00,
    771.00,
    829.00,
    971.00,
   1029.00,
   1271.00,
   1329.00,
   1571.00,
   1629.00,
   1871.00,
   2129.00,
   2371.00,
   2629.00])

  kzi = ZI.shape[0]
#
# Follow ncoda_archv_inc.f to define the mid-depths points in NCODA
# Note that in the ncoda code zz is NCODA z-level mid-depths points
# and zi is array with z-level interface depths
  ZM = np.zeros((kzi-1))
  for kk in range(1,kzi-1):
    ZM[kk] = 0.5*(ZI[kk-1]+ZI[kk])

  kzm = ZM.shape[0]

  if zneg:
    ZI = -ZI
    ZM = - ZM

  return ZI, kzi, ZM, kzm

def find_indx_lonlat(x0, y0, X0, Y0, dlt_err=7000., fatal_err=True):
  """
  Given pnt coordinates (x0,y0) and 2D lon/lat arrays (X0, Y0):
  Find closest grid point (ii0,jj0) to lon/lat coordinate
  For W-E sections - use find_indx_2Dsect 

  dlt_err - max difference (degrees) between x0,y0 and found point on X0,Y0 before warning/error
  fatal_error - fatal, otherwise - warning with indices < 0 returned 

  Output: ii0, jj0 - indices closest to x0,y0 (lon, lat)
  """
  import mod_misc1 as mmisc
  
  if abs(y0) > 90.:
    raise Exception(f'latitude y0 >/< 90N/S, check y0 is lat: {y0}')

  # Normalize x0 or shift longitudes wrt to x0:
  #x0 = (x0 + 360) % 360
  #X0 = (X0 + 360) % 360
  
  # Shift longitudes wrt to x0
  XX = ((X0 - x0 + 180) % 360) - 180 + x0
  #XX = X0.copy()
  YY = Y0.copy()

  #if np.max(XX) > 180.:
  #  XX = np.where(XX > 180., XX-360., XX)

  dmm = mmisc.dist_sphcrd(y0,x0,YY,XX)   # distance in meters
  jj0, ii0 = np.unravel_index(np.argmin(dmm), dmm.shape)
  dmm_min = np.min(dmm)

# Sanity check:
#  print(f'jj0={jj0} ii0={ii0} x0={x0} y0={y0}')
  derr = mmisc.dist_sphcrd(y0,x0,YY[jj0,ii0],XX[jj0,ii0])   # error in meters
  #derr = np.sqrt((XX[jj0,ii0]-x0)**2 + (YY[jj0,ii0]-y0)**2)
  if derr >= dlt_err:
    ii0 = -ii0
    jj0 = -jj0
    print(f'WARNING: Found index is far away from given from lon/lat: Error = {derr*1e-3} km')
    print(f'WARNING: Check lon0={x0} lat0={y0} matched to lon={XX[jj0,ii0]} lat={YY[jj0,ii0]}')
    if fatal_err:
      raise Exception(f'ERROR Finding index from lon/lat')

  return ii0, jj0

def find_indx_2Dsect(x0,y0,X0,Y0,xsct="none", dlt_err=1, fatal_err=True):
  """
  Given pnt coordinates (x0,y0) and 2D lon/lat arrays (X0, Y0):
  Find closest grid point (ii0,jj0) to lon/lat coordinate
  For W-E sections, provide xsct name 
  then indices are given relative to the section 1st index

  dlt_err - max difference (degrees) between x0,y0 and found point on X0,Y0 before warning/error
  fatal_error - fatal, otherwise - warning with indices < 0 returned 

  Output: ii0, jj0 - indices closest to x0,y0 (lon, lat)
          ip0, jp0 - indices relative to 1st pnt in the section
  """
  import mod_misc1 as mmisc
  if x0 > 180.:
    x0 = x0-360.

  if abs(y0) > 90.:
    raise Exception(f'latitude y0 >/< 90N/S, check y0 is lat: {y0}')

  XX = X0.copy()
  YY = Y0.copy()

  if np.max(XX) > 180.:
    XX = np.where(XX > 180., XX-360., XX)

  dmm = np.sqrt((XX-x0)**2+(YY-y0)**2)
  jj0, ii0 = np.unravel_index(np.argmin(dmm), dmm.shape)
  dmm_min = np.min(dmm)
#
# Check for singularity along 180/-180 longitude
  IDM = XX.shape[1]
  if ii0 > 0 and ii0 < IDM-1:
    xm1 = XX[jj0,ii0-1]
    xp1 = XX[jj0,ii0+1]
    if abs(xm1-x0) > 180. or abs(xp1-x0)>180.:
      XX = np.where(XX < 0., XX+360., XX)
      if x0 < 0:
        x0 = x0+360.
    dmm = np.sqrt((XX-x0)**2+(YY-y0)**2)

    if np.min(dmm) < dmm_min:
      jj0, ii0 = np.where(dmm == np.min(dmm)) # global indices
      jj0 = jj0[0]
      ii0 = ii0[0]

  if ii0 == IDM-1:
    xm1 = XX[jj0,ii0]
    xp1 = XX[jj0,0]
    if abs(xp1-x0) > 180.:
      XX = np.where(XX < 0., XX+360., XX)
      if x0 < 0:
        x0 = x0+360.
    dmm = np.sqrt((XX-x0)**2+(YY-y0)**2)

    if np.min(dmm) < dmm_min:
      jj0, ii0 = np.where(dmm == np.min(dmm)) # global indices
      jj0 = jj0[0]
      ii0 = ii0[0]

# Sanity check:
#  print(f'jj0={jj0} ii0={ii0} x0={x0} y0={y0}')
  derr = mmisc.dist_sphcrd(y0,x0,YY[jj0,ii0],XX[jj0,ii0])   # error in meters
  #derr = np.sqrt((XX[jj0,ii0]-x0)**2 + (YY[jj0,ii0]-y0)**2)
  if derr >= dlt_err:
    ii0 = -ii0
    jj0 = -jj0
    print(f'WARNING: Found index is far away from given from lon/lat: Error = {derr}')
    print(f'WARNING: Check lon0={x0} lat0={y0} matched to lon={XX[jj0,ii0]} lat={YY[jj0,ii0]}')
    if fatal_err:
      raise Exception(f'ERROR Finding index from lon/lat')
#
# If section provided, indices for this section:
  ip0 = -1
  jp0 = -1

  if not xsct=="none":
    SCT = Jsections()
    aa  = SCT.get(xsct)
    i1  = aa[0]
    i2  = aa[1]
    j1  = aa[2]
    ip0 = ii0-i1
    jp0 = jj0-j1

    return ii0, jj0, ip0, jp0
  else:
    return ii0, jj0


def interp_indx_lonlat(x0,y0,LON,LAT):
  """
    Map x0,y0 (lon,lat) into index space, i.e.
    find "exact" index (float number) by interpolation
  """
  import mod_bilinear as mblinr
#
# Find 4 nodes around x0,y0
# arranging counter-clockwise 
  i1, j1 = find_indx_lonlat(x0,y0,LON,LAT)
  x1 = LON[j1,i1]
  y1 = LAT[j1,i1]

  JDM = LON.shape[0]
  IDM = LON.shape[1]
 
  if x1 <= x0:
    iv1 = i1
    iv2 = i1+1
  else:
    iv1 = i1-1
    iv2 = i1

  if y1 <= y0:
    jv1 = j1
    jv2 = j1+1
  else:
    jv1 = j1-1
    jv2 = j1

  iht1 = iv1  # for interpolation keep original values
  iht2 = iv2
  jht1 = jv1
  jht2 = jv2
  if iv1 < 0:
    iv1 = IDM-1
  if jv1 < 0:     # impossible case
    jv1 = 0
  if iv2 > IDM-1:
    iv2 = 0
  if jv2 > JDM-1:
    jv2 = JDM-1
#
# Form 1D arrays of X and Y coordinates of vertices
# of a grid cell including the interp point
  XX = np.array([LON[jv1,iv1],LON[jv1,iv2],LON[jv2,iv2],LON[jv2,iv1]])
  YY = np.array([LAT[jv1,iv1],LAT[jv1,iv2],LAT[jv2,iv2],LAT[jv2,iv1]])
  
# Singularity:
  if (np.max(XX)-np.min(XX)) > 180.:
    II =  np.where(XX<0)
    XX[II] = XX[II]+360.
    if x0 < 0.:
      x0 = x0+360.

# Vector of indices:
  IX = np.array([iht1,iht2,iht2,iht1])
  JX = np.array([jht1,jht1,jht2,jht2])

# Map X,Y ---> Xhat, Yhat on reference quadrialteral 
# i.e. map WOA grid coordinate to a reference quadrilateral 
# to do bilinear interpolation 
  xht, yht = mblinr.map_x2xhat(XX,YY,x0,y0)

# Perform interpolation on reference rectangle, that is 
# similar to interp on actual rectangle
  Iint = mblinr.blinrint(xht,yht,IX)
  Jint = mblinr.blinrint(xht,yht,JX)

# Check negative values:
  if Iint < 0:
    Iint = IDM+Iint
  elif Iint > IDM-1:
    Iint = Iint-(IDM-1)

  return Iint, Jint


def anls_lrthkn(ZZ,lrintrf=True):
  """
  Analyze layer interface depths for large jumps
  ZZ - default, interface depths lrintrf
       if not - layer thickness (m)
  """
  if not lrintrf:
    dZZ = ZZ.copy()
  else:
    dZZ = np.abs(np.diff(ZZ, axis=0))

  kdm = dZZ.shape[0]
  jdm = dZZ.shape[1]
  idm = dZZ.shape[2]
 
  grdZ = np.zeros((jdm,idm))
  for k in range(kdm):
    aa = np.squeeze(dZZ[k,:,:])
#    aa = np.where(np.isnan(aa),0,aa)
    dzdx = np.zeros((jdm,idm))
    dzdy = np.zeros((jdm,idm))
    for i in range(1,idm):
      dzdx[:,i] = np.abs(aa[:,i-1]-aa[:,i])
    for j in range(1,jdm):
      dzdy[j,:] = np.abs(aa[j-1,:]-aa[j,:])

    dmm   = dzdx-dzdy
    gradz = np.where(dmm>0,dzdx,dzdy)
    dmm   = gradz-grdZ
    grdZ  = np.where(dmm>0,gradz,grdZ)

    print(' Grad lr thkn: k={0} max={1:8.2f} m'.format(k,np.nanmax(grdZ)))
  
  return grdZ

def find_max3d(A3d,dZZ):
  """
  Find max values over all layers for 3D array
  avoid zero thickness layers
  """
  kdm = A3d.shape[0]
  jdm = A3d.shape[1]
  idm = A3d.shape[2]
 
  maxA = np.zeros((jdm,idm))-2.e20
  for k in range(kdm):
    aa    = np.squeeze(A3d[k,:,:])
    dz    = np.squeeze(dZZ[k,:,:])
    dmm   = aa-maxA
    dmm   = np.where(dz<1.e-6,-1.e3,dmm)
    maxA  = np.where(dmm>0.,aa,maxA) 
    print('Max value: k={0} max={1:10.4f}'.format(k,np.nanmax(aa)))
  
  return maxA

def find_min3d(A3d,dZZ):
  """
  Find min values over all layers for 3D array
  avoid zero thickness layers
  """
  kdm = A3d.shape[0]
  jdm = A3d.shape[1]
  idm = A3d.shape[2]
 
  minA = np.zeros((jdm,idm))+2.e20
  for k in range(kdm):
    aa    = np.squeeze(A3d[k,:,:])
    dz    = np.squeeze(dZZ[k,:,:])
    dmm   = minA-aa
    dmm   = np.where(dz<1.e-6,-1.e3,dmm)  # zero-thickness layers
    minA  = np.where(dmm>0.,aa,minA) 
    print('Min value: k={0} min={1:10.4f}'.format(k,np.nanmin(aa)))
  
  return minA


def clrmp_lmask(nclrs=2,clr_land=[0.3,0.3,0.3]):
  """
    Create colormap for land mask with 2 colors
  """

  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  r, g, b = clr_land[0:3]
  clrs   = cm.get_cmap('GnBu_r',nclrs)
  newclr = clrs(range(nclrs))
  newclr[0,:] = [1, 1, 1, 1]
  newclr[1,:] = [r, g, b, 1]

  newcmp  = ListedColormap(newclr)

  return newcmp


def interp_1Dhycom(ZS,SS,ZZh,ZMh):
  """
    Interp SS profile at ZS depths onto
    HYCOM mid-point depths ZMh
    depths are negative
  """
  import mod_interp1D as mintrp

  if np.min(ZS) > 0.:
    ZS = -ZS
  if np.min(ZMh) > 0.:
    ZMh = -ZMh
  if np.min(ZZh) > 0.:
    ZZh = -ZZh

# Add missing near-surf values in ARGO
  if abs(ZS[0]) > 1.e-3:
    ZS = np.insert(ZS,0,0.)
    SS = np.insert(SS,0,SS[0])

  kzm = ZMh.shape[0]
  kzz = ZZh.shape[0]
  ksz = ZS.shape[0]
  iSS = np.zeros((kzm))*np.nan
  for kk in range(kzm):
    z1 = ZZh[kk]
    z2 = ZZh[kk+1]
    z0 = ZMh[kk]
    if np.isnan(z0):
      break

    if min(ZS) > z2:   # Argo last point is too shallow
      break
#
# Interpolate SS onto upper layers interface
    isz1 = max(np.where(ZS >= z1)[0])
    isz2 = min(np.where(ZS <= z2)[0])

    if isz1 == isz2:
      iSS[kk] = SS[isz1]
      continue

# Prepare Argo T/S and depths for HYCOM layer
    zz  = ZS[isz1:isz2+1]
    xx  = SS[isz1:isz2+1]
    xx1 = mintrp.pcws_lagr1(zz,xx,z1)
    xx2 = mintrp.pcws_lagr1(zz,xx,z2)
# 
# Integrate Argo over the HYCOM layer using trapezoidal rule:
    zsgm  = zz[1:-1]
    zsgm  = np.insert(zsgm,0,z1)
    zsgm  = np.append(zsgm,z2)
    dzsgm = abs(np.diff(zsgm))
    xsgm  = xx[1:-1]
    xsgm  = np.insert(xsgm,0,xx1)
    xsgm  = np.append(xsgm,xx2)
    lsgm  = np.shape(dzsgm)[0]
    asum  = 0.
    for ll in range(lsgm):
      Ia = 0.5*dzsgm[ll]*(xsgm[ll]+xsgm[ll+1])
      asum += Ia
    dZ = np.sum(dzsgm)  # should be = z1-z2 - HYCOM layer
    iSS[kk] = asum/dZ


#    zz0 = np.arange(zz[0],zz[-1],-0.1)
#    A = []
#    for ll in range(zz0.shape[0]):
#      x0 = mintrp.pcws_lagr1(zz,xx,zz0[ll])
#      A.append(x0)
#    A = np.array(A)
#    plt.plot(A,zz0)

  return iSS

def err_stat(S1,S2):
  """
    Compute l2 norm, RMSE, inf-norm
    for 2 vectors S1, S2
  """
  r1  = (S1-S2)**2
  inn = np.where(~np.isnan(r1))[0][-1]

  Sl2   = np.sqrt(np.dot(r1[:inn+1],r1[:inn+1]))
  Srmse = np.sqrt(np.dot(r1[:inn+1],r1[:inn+1])/inn)
  Sinf  = np.nanmax(abs(S1-S2))

  return Sl2, Srmse, Sinf

def prof_limits(Th,ZMh,nT=5,z0=-6000.,ndec=0):
  """
    Define min/max limits for 1D profile Th, ZM - depths
    over depth from z0:surface
    nT - approximate number of x ticks
    ndec - # of decimals for rounding when search for min/max limits
  """
  iz  = max(np.where(ZMh >= z0)[0])+1
  dmm = Th[:iz+1]
  dmx = np.max(dmm)
  dmn = np.min(dmm)
  dT  = int(abs(dmx-dmn))
  if dT == 0:
    dT = 1

  cff = 10.**ndec
  tlim1 = np.floor(dmn)
  tlim2 = np.ceil(dmx)
  if ndec > 0:
    tlim1 = 1./cff*(np.floor(dmn*cff))
    tlim2 = 1./cff*(np.ceil(dmx*cff))

  if tlim1 == tlim2:
    cff2 = 10.*cff
    tlim1 = 1/cff2*(np.floor(dmn*cff2))
    tlim2 = 1/cff2*(np.ceil(dmx*cff2))

  dltT = abs(tlim2-tlim1)/nT

  return tlim1, tlim2, dltT

def axis_limits2(T1,T2,cff=1):
  """
    Define axis limits for 2 1D arrays
    cff - rounding precision, 1- intigers, 0.1, ...etc
  """
  t1 = min([np.min(T1),np.min(T2)])
  t2 = max([np.max(T1),np.max(T2)])

  tlim1 = cff*np.floor(t1/cff)
  tlim2 = cff*np.ceil(t2/cff)

  return tlim1, tlim2
  
def pool_data_zbins(T1, T2, ZM, Zbins):
  """
    Pool data from 2 profiles by depth bins
    for scatter plots
    Profiles are at the same depths 
    ZM < 0
  """

  class tbinned():
    kind = 'Binned data'

    def __init__(self):
      self.T1b = []
      self.T2b = []

    def add_data(self,T1b,T2b):
      self.T1b.append(T1b)       # binned data from prof 1
      self.T2b.append(T2b)       # binned data from prof 2

  TB = tbinned()
  nbins = Zbins.shape[0]-1
 
  zmin = np.nanmin(ZM) 
  for kk in range(nbins):
    z1 = Zbins[kk]
    z2 = Zbins[kk+1]
    if z2 > z1:
      z1 = Zbins[kk+1]
      z2 = Zbins[kk]
#
# Profile is too shallow for depth bin
    if zmin > z1:
      dmm1 = np.nan
      dmm2 = np.nan
    else:
      iz1 = min(np.where(ZM <= z1)[0])
      iz2 = max(np.where(ZM > z2)[0])

      dmm1 = T1[iz1:iz2+1]
      dmm2 = T2[iz1:iz2+1]

    TB.add_data(dmm1,dmm2)

  return TB

class tserr():
  kind = 'TS Argo error stat'

  def __init__(self):
    self.recn  = []
    self.numb  = []
    self.lon   = []
    self.lat   = []
    self.ptime = []
    self.rtime = []
    self.Tstd  = []
    self.Sstd  = []
    self.SL2   = []
    self.TL2   = []
    self.SRMSE = []
    self.TRMSE = []
    self.SINF  = []
    self.TINF  = []
    self.SS01  = []  # pooled S values in the top level 01
    self.TT01  = []
    self.SS02  = []
    self.TT02  = []
    self.Sqcrj = []  # indices of rejected S profiles
    self.Tqcrj = []  # omdoces pf rejected T profiles

  def add_data(self, recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
               SL2, TL2, SRMSE, TRMSE, SINF, TINF, SS01, TT01, \
               SS02, TT02, Sqcrj, Tqcrj):
    self.recn.append(recn)    # record #
    self.numb.append(numb)    # Argo #
    self.lon.append(lon)
    self.lat.append(lat)
    self.ptime.append(dnmb)
    self.rtime.append(rnmb)
    self.Tstd.append(Tstd)
    self.Sstd.append(Sstd)
    self.SL2.append(SL2)
    self.TL2.append(TL2)
    self.SRMSE.append(SRMSE)
    self.TRMSE.append(TRMSE)
    self.SINF.append(SINF)
    self.TINF.append(TINF)
    self.SS01.append(SS01)     #values in the top level 01
    self.TT01.append(TT01)
    self.SS02.append(SS02)
    self.TT02.append(TT02)
    self.Sqcrj.append(Sqcrj) #S rejected profiles
    self.Tqcrj.append(Tqcrj)



class tsarray():
  def __init__(self,recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
               SL2, TL2, SRMSE, TRMSE, SINF, TINF, SS01, TT01, \
               SS02, TT02, Sqcrj, Tqcrj):
    self.recn  = recn
    self.numb  = numb
    self.lon   = lon
    self.lat   = lat
    self.ptime = dnmb
    self.rtime = rnmb
    self.Tstd  = Tstd
    self.Sstd  = Sstd
    self.SL2   = SL2
    self.TL2   = TL2
    self.SRMSE = SRMSE
    self.TRMSE = TRMSE
    self.SINF  = SINF
    self.TINF  = TINF
    self.SS01  = SS01
    self.TT01  = TT01
    self.SS02  = SS02
    self.TT02  = TT02
    self.Sqcrj = Sqcrj
    self.Tqcrj = Tqcrj

  def add_array(self, recn, numb, lon, lat, dnmb, rnmb, Tstd, Sstd, \
               SL2, TL2, SRMSE, TRMSE, SINF, TINF, SS01, TT01, \
               SS02, TT02, Sqcrj, Tqcrj):
    self.recn  = np.append(self.recn,recn, axis=0)
    self.numb  = np.append(self.numb,numb, axis=0)
    self.lon   = np.append(self.lon,lon, axis=0)
    self.lat   = np.append(self.lat,lat, axis=0)
    self.ptime = np.append(self.ptime,dnmb, axis=0)
    self.rtime = np.append(self.rtime,rnmb, axis=0)
    self.Tstd  = np.append(self.Tstd,Tstd, axis=0)
    self.SL2   = np.append(self.SL2,SL2, axis=0)
    self.TL2   = np.append(self.TL2,TL2, axis=0)
    self.SRMSE = np.append(self.SRMSE,SRMSE, axis=0)
    self.TRMSE = np.append(self.TRMSE,TRMSE, axis=0)
    self.SINF  = np.append(self.SINF,SINF, axis=0)
    self.TINF  = np.append(self.TINF,TINF, axis=0)
    self.SS01  = np.append(self.SS01,SS01, axis=0)
    self.TT01  = np.append(self.TT01,TT01, axis=0)
    self.SS02  = np.append(self.SS02,SS02, axis=0)
    self.TT02  = np.append(self.TT02,TT02, axis=0)
    self.Sqcrj = np.append(self.Sqcrj,Sqcrj, axis=0)
    self.Tqcrj = np.append(self.Tqcrj,Tqcrj, axis=0)


class med_iqr():
  def __init__(self, med, iqrbt, iqrup, p10, p90, vmin, vmax):
    self.median = np.array([med], dtype=float)
    self.iqrbt  = np.array([iqrbt], dtype=float)
    self.iqrup  = np.array([iqrup], dtype=float)
    self.p10    = np.array([p10], dtype=float)
    self.p90    = np.array([p90], dtype=float)
    self.vmin   = np.array([vmin], dtype=float)
    self.vmax   = np.array([vmax], dtype=float)

  def add_data(self, med, iqrbt, iqrup, p10, p90, vmin, vmax):
    self.median = np.append(self.median, med)
    self.iqrbt  = np.append(self.iqrbt, iqrbt)
    self.iqrup  = np.append(self.iqrup, iqrup)
    self.p10    = np.append(self.p10, p10)
    self.p90    = np.append(self.p90, p90)
    self.vmin   = np.append(self.vmin, vmin)
    self.vmax   = np.append(self.vmax, vmax)
 

class rmmx():
  def __init__(self,KDM):
    self.rmin   = np.zeros((KDM,1), dtype=float)
    self.rmax   = np.zeros((KDM,1), dtype=float)
    self.dh0max = np.zeros((KDM,1), dtype=float)
    self.dh0min = np.zeros((KDM,1), dtype=float)
    self.dhmax  = np.zeros((KDM,1), dtype=float)
    self.dhmin  = np.zeros((KDM,1), dtype=float)
    self.imin   = np.zeros((KDM,1), dtype=int)
    self.jmin   = np.zeros((KDM,1), dtype=int)
    self.imax   = np.zeros((KDM,1), dtype=int)
    self.jmax   = np.zeros((KDM,1), dtype=int)
    self.time   = np.zeros((1), dtype=int)

  def add_data(self, rmin, rmax, imin, jmin, imax, jmax, \
                dh0min, dh0max, dhmin, dhmax, dnmb):
    self.rmin   = np.append(self.rmin, rmin, axis=1)
    self.rmax   = np.append(self.rmax, rmax, axis=1)
    self.dh0max = np.append(self.dh0max, dh0max, axis=1)
    self.dh0min = np.append(self.dh0min, dh0min, axis=1)
    self.dhmax  = np.append(self.dhmax, dhmax, axis=1)
    self.dhmin  = np.append(self.dhmin, dhmin, axis=1)
    self.imin   = np.append(self.imin, imin, axis=1)
    self.jmin   = np.append(self.jmin, jmin, axis=1)
    self.imax   = np.append(self.imax, imax, axis=1)
    self.jmax   = np.append(self.jmax, jmax, axis=1)
    self.time   = np.append(self.time, dnmb)



def unpack_TSbins_Iprf(dmm,Ipndx,ibin):
  """
    Pool profiles for selected Ipndx profiles 
    for depth bin ibin
    into 1 array for error stat or plotting
    dmm is an array nrec x nbins
    nrec = # of profiles
    nbins = # of depth bins 
    dmm - same type of obs/model profiles (Argo, RTOFS, ...)
  """
  nrec  = dmm.shape[0]
  nbins = dmm.shape[1]
  tbb   = dmm[:,ibin]  
  nslct = Ipndx.shape[0]

  for kk in range(nslct):
    ii = Ipndx[kk]
    ai = tbb[ii]
    if isinstance(ai, np.ndarray):  # missing depths in profiles, ai=np.nan
      if kk == 0:
        Prf = ai.copy()
      else:
        Prf = np.append(Prf,ai, axis=0)
    else:
      continue

  return Prf


def unpack_TSbins(dmm,ibin):
  """
    Pool all profiles for depth bin ibin
    into 1 array for plotting
    dmm is an array nrec x nbins
    nrec = # of profiles
    nbins = # of depth bins 
    dmm - same type of obs/model profiles (Argo, RTOFS, ...)
  """
  nrec  = dmm.shape[0]
  nbins = dmm.shape[1]
  tbb   = dmm[:,ibin]  

  for ii in range(nrec):
    ai = tbb[ii]
    if isinstance(ai, np.ndarray):  # missing depths in profiles
      if ii == 0:
        Prf = ai.copy()
      else:
        Prf = np.append(Prf,ai, axis=0)
    else:
      continue

  return Prf


def find_rejected_TSbinned(Irj,dmm,ibin):
  """
    Find rejected TS values of binned data
    Irj - indices of rejected T/S profiles in
          unordered arrays packed in structures
    dmm - group of T/S data in depth bins,  ndarray
          dmm is an ndarray N x nbins
          N arrays (N is the total # of profiles)
          nbins = # of depth bins
          each array is a group of T /S data in this depth bin
    Need to combine all T/S of rejected profiles into 1 array/list
  """
  nrec  = dmm.shape[0]
  nbins = dmm.shape[1]
  tbb   = dmm[:,ibin]
  nrjct = Irj.shape[0]
  Trjct = []
  for ii in range(nrjct):
    ixx = Irj[ii]
    Trjct.append(tbb[ixx])

  return Trjct

def find_notrejected_TSbinned(Irj,dmm,ibin):
  nrec  = dmm.shape[0]
  nbins = dmm.shape[1]
  tbb   = dmm[:,ibin]
  nrjct = Irj.shape[0]
  Tnrjct = []
  for ii in range(nrec):
    if not np.any(Irj == ii):
      Tnrjct.append(tbb[ii])

  return Tnrjct
  

def find_val_TSbinned(dmm,vmin,vmax,ibin,Irj,Fnrj=True):
  """
  Find value within [vmin,vmax] for not rejected (Fnrj=True) 
  or rejected (False) values
  dmm is either T or S from TSERR structured object
  for N profiles x N depth bins
  e.g., (470, 2)
  includes both rejected and not rejected profiles
  ibin - depth bin (0, 1, ...)
  Irj - indices of rejected profiles that should be discarded in the search
  """
  nrec  = dmm.shape[0]
  nbins = dmm.shape[1]
  tbb   = dmm[:,ibin]

#  Tnrj = find_notrejected_TSbinned(Irj,dmm,ibin)
#  nnrj = len(Tnrj)
  Indx = []
  for ii in range(nrec):
    if Fnrj:
      if np.any(Irj == ii):
        continue

    else:
      if not np.any(Irj == ii):
        continue

    aa = tbb[ii]
    imm = np.where( (aa >= vmin) & (aa <= vmax) )[0]
    if len(imm) > 0:
      Indx.append(ii)

  return Indx

def plot_boxplot(axs,ST,ctl='boxplot',Xlbls=[]):
  """
  Plot boxplot on axes axs
  all statistics are in the ST structure array
  """
  clrbx = [0.,0.4,0.9]
  clrmd = [0.8,0.1,0]
  clrln = [0.4,0.4,0.4]
  dbx   = 0.2   # half of box width
  dtg   = 0.05  # hash tag
  Nbins = len(ST.median)+1

  ymin = 1.e6
  ymax = -1.e6
  for kk in range(Nbins-1):
    mdn  = ST.median[kk]
    iq1  = ST.iqrbt[kk]
    iq2  = ST.iqrup[kk]
    p10  = ST.p10[kk]
    p90  = ST.p90[kk]
    vmin = ST.vmin[kk]
    vmax = ST.vmax[kk]

    x0  = kk+1
    xb1 = x0-dbx
    xb2 = x0+dbx
    xt1 = x0-dtg
    xt2 = x0+dtg

    axs.plot([xb1,xb2],[mdn,mdn],color=clrmd)
    axs.plot([xb1,xb2],[iq1,iq1],color=clrbx)
    axs.plot([xb1,xb2],[iq2,iq2],color=clrbx)
    axs.plot([xb1,xb1],[iq1,iq2],color=clrbx)
    axs.plot([xb2,xb2],[iq1,iq2],color=clrbx)
    axs.plot([x0,x0],[p10,iq1],color=clrln)
    axs.plot([x0,x0],[iq2,p90],color=clrln)
    axs.plot([xt1,xt2],[p10,p10],color=clrln)
    axs.plot([xt1,xt2],[p90,p90],color=clrln)
  #  axs.plot([0,Nbins-1],[0,0],'--',color=[0.5

    ymin = min([ymin,p10])
    ymax = max([ymax,p90])


  dyy = abs(ymax-ymin)
  yl1 = 0.1*(np.floor((ymin-0.05*dyy)*10))
  yl2 = 0.1*(np.ceil((ymax+0.05*dyy)*10))
  if yl1 < 0. and yl2 > 0.:
    ylim = max([abs(yl1),abs(yl2)])
    yl1 = -ylim
    yl2 = ylim

  dyl  = abs(yl2-yl1)
  ntck = 10
  dtck = dyl/ntck
  TCK = np.array([0.001,0.025,0.05,0.1,0.25,0.5,1.,2.5,5.,10.,25.,50.,100.])
  tmm = abs(TCK-dtck)
  it  = np.where(tmm == np.min(tmm))[0]
  dtck = TCK[it]

  axs.set_xticks(np.arange(1,Nbins))
  axs.set_yticks(np.arange(np.floor(yl1),np.ceil(yl2),dtck))
  axs.set_ylim(yl1,yl2)
  if len(Xlbls) > 0:
    axs.set_xticklabels(Xlbls)
  axs.grid(True)
  axs.set_title(ctl)
    
  return axs

def plot_boxplot_v2(ax1, AA, sttl='boxplots', CLRS=[], XX=[], clrln=[0,0,0], \
                    xlbl='x axis', ylbl='y axis', lgnd_names=[], btx=[]):
  """
    Plot boxplot diagr for 3D AA array [Ngroups x stat x data]
    stat = median, 25th, 75th percentiles, min, max values for whiskies
           5 values
    Can be 2D (if N = 1) array, will be expanded to 3D

    Note correct shape of AA, an option of different order of the axes
    can be added, i.e. allow stat x data x Ngroups
    
    Provide lgnd_names for legend
    lgnd_names = list of strings, len(lgnd_names) = Ngrps

  """
  import mod_colormaps as mclrmps
  from matplotlib.patches import Polygon

  ndim = len(AA.shape)
  if ndim == 2: 
    AA = np.expand_dims(AA, axis=2)
  kdm, jdm, idm = AA.shape

  Ngrps = kdm 

  if len(CLRS) == 0:
    CLRS = mclrmps.colormap_discrete(cmp_obj=False)
    if Ngrps > 10:
      CLRS = mclrmps.create_colormap(CLRS, Ngrps, cmp_obj=False)

# Define box width and distance between boxes:
  if len(XX) == 0:
    XX = np.arange(idm)
# Add some space between boxplots along XX axis 
# for visual separation: d1box = 1/Ngrps - all boxplots have
# sames separation, 0.8/Ngrps - boxplots are sqiuzzed within XX intervals
  dend    = 0.2
  d1box   = (1.-2*dend)/Ngrps  # X space for 1 boxplot with spaces with spaces at the ends
  dspace  = 0.1*d1box  # half-spacer between boxes
  dbx_hlf = 0.5*(d1box - 2*dspace)  # half box width
  dtg   = 0.05  # hash tag
   
#  plt.clf()
#  ax1 = plt.axes([0.1, 0.5, 0.8, 0.4])
  for imo in range(idm):
    for igrp in range(Ngrps):
      xx0 = XX[imo] +dend + igrp*(2*dspace + 2*dbx_hlf) + dspace + dbx_hlf
      stat = AA[igrp,:,imo].squeeze()
      if len(stat) >= 5:
        mdn, lprc, uprc, tmin, tmax = stat[:5]

      verts = [(xx0-dbx_hlf, lprc),(xx0-dbx_hlf,uprc),(xx0+dbx_hlf,uprc),\
               (xx0+dbx_hlf,lprc)]
      fclr = CLRS[igrp]
      poly = Polygon(verts, facecolor=fclr)
      ax1.add_patch(poly)

      # Draw a box
      ax1.plot([xx0-dbx_hlf, xx0-dbx_hlf],[lprc,uprc],'-',linewidth=1, color=clrln)
      ax1.plot([xx0+dbx_hlf, xx0+dbx_hlf],[lprc,uprc],'-',linewidth=1, color=clrln)
      ax1.plot([xx0-dbx_hlf, xx0+dbx_hlf],[uprc,uprc],'-',linewidth=1, color=clrln)
      ax1.plot([xx0-dbx_hlf, xx0+dbx_hlf],[lprc,lprc],'-',linewidth=1, color=clrln)

      # Draw whiskies showing the range of the data
      ax1.plot([xx0-dbx_hlf,xx0+dbx_hlf],[mdn,mdn],linewidth=2, color=clrln)
      ax1.plot([xx0,xx0],[uprc,tmax],'-',linewidth=2, color=clrln)
      ax1.plot([xx0,xx0],[tmin,lprc],'-',linewidth=2, color=clrln)
      ax1.plot([xx0-dtg,xx0+dtg],[tmax,tmax],'-',linewidth=2, color=clrln)
      ax1.plot([xx0-dtg,xx0+dtg],[tmin,tmin],'-',linewidth=2, color=clrln)

  ax1.set_xticks(XX)
  ax1.set_xlim([XX[0],XX[-1]+1])
  ax1.grid('on')
  ax1.set_xlabel(xlbl)
  ax1.set_ylabel(ylbl)
  ax1.set_title(sttl)

  if len(lgnd_names) > 0:
    assert len(lgnd_names) == Ngrps, f'Number of legend names should be {Ngrps}'
    ax2 = plt.axes([0.6, 0.1, 0.45, 0.12])
    dltY = 0.015
    for igrp in range(Ngrps):
      y0 = igrp*dltY + dltY
      x0 = 0.1
      fclr = CLRS[igrp]
      grp_name = lgnd_names[igrp]
      ax2.plot(x0, y0, marker='s', ms=15, markerfacecolor=fclr, mec='none')
      ax2.text(x0+0.1, y0, grp_name)

    ax2.set_ylim([0, y0+dltY])
    ax2.set_xlim([0, 0.8])
    ax2.axis('off')

  if len(btx) > 0:
    bottom_text(btx, pos=[0.05, 0.08])

  plt.sca(ax1)

  return ax1

def plot_boxplot2D(ax1, AA, sttl='boxplots', CLRS=[], XX=None, clrln=[0,0,0], \
                    xlbl='x axis', ylbl='y axis', grp_names=[], btx=[]):
  """
    Plot boxplot diagr for 2D AA array [Ngroups x stat]
    stat = median, 25th, 75th percentiles, min, max values for whiskies
           5 values

    Provide grp_names for X labeling
    grp_names = list of strings, len(grp_names) = Ngrps

  """
  import mod_colormaps as mclrmps
  from matplotlib.patches import Polygon

  ndim = len(AA.shape)
  assert ndim == 2, 'Stat array should be 2D: groups x stat data'
  kdm, jdm = AA.shape
  assert jdm >= 5, 'The number of statistics should be >=5: median, 25,75 prc, min, max'

  Ngrps = kdm 

  if len(CLRS) == 0:
    CLRS = np.zeros((Ngrps,3))
    for kk in range(Ngrps):
      CLRS[kk,:] = np.array([0.6,0.6,0.6])

# Define box width and distance between boxes:
  if XX is None or len(XX) == 0:
    XX = np.arange(1,Ngrps+1)

# Add some space between boxplots along XX axis 
# for visual separation: d1box = 1/Ngrps - all boxplots have
# sames separation, 0.8/Ngrps - boxplots are sqiuzzed within XX intervals
  d1box   = 0.6     # box widths
  dspace  = 0.1*d1box  # half-spacer between boxes
  dbx_hlf = 0.5*(d1box - 2*dspace)  # half box width
  dtg   = 0.05  # hash tag
   
  for igrp in range(Ngrps):
    xx0 = XX[igrp]
    stat = AA[igrp,:].squeeze()
    if len(stat) >= 5:
      mdn, lprc, uprc, tmin, tmax = stat[:5]

    verts = [(xx0-dbx_hlf, lprc),(xx0-dbx_hlf,uprc),(xx0+dbx_hlf,uprc),\
             (xx0+dbx_hlf,lprc)]
    fclr = CLRS[igrp]
    poly = Polygon(verts, facecolor=fclr, zorder=2)
    ax1.add_patch(poly)

    # Draw a box
    ax1.plot([xx0-dbx_hlf, xx0-dbx_hlf],[lprc,uprc],'-',linewidth=1, color=clrln)
    ax1.plot([xx0+dbx_hlf, xx0+dbx_hlf],[lprc,uprc],'-',linewidth=1, color=clrln)
    ax1.plot([xx0-dbx_hlf, xx0+dbx_hlf],[uprc,uprc],'-',linewidth=1, color=clrln)
    ax1.plot([xx0-dbx_hlf, xx0+dbx_hlf],[lprc,lprc],'-',linewidth=1, color=clrln)

    # Draw whiskies showing the range of the data
    ax1.plot([xx0-dbx_hlf,xx0+dbx_hlf],[mdn,mdn],linewidth=2, color=clrln)
    ax1.plot([xx0,xx0],[uprc,tmax],'-',linewidth=2, color=clrln)
    ax1.plot([xx0,xx0],[tmin,lprc],'-',linewidth=2, color=clrln)
    ax1.plot([xx0-dtg,xx0+dtg],[tmax,tmax],'-',linewidth=2, color=clrln)
    ax1.plot([xx0-dtg,xx0+dtg],[tmin,tmin],'-',linewidth=2, color=clrln)

  if len(XX)>1:
    dX = XX[1]-XX[0]
  else:
    dX = 0.8

  ax1.set_xticks(XX)
  ax1.set_xlim([XX[0]-dX,XX[-1]+dX])
  ax1.grid('on')
  ax1.set_xlabel(xlbl)
  ax1.set_ylabel(ylbl)
  ax1.set_title(sttl)

  if len(grp_names) == len(XX):
    ax1.set_xticklabels(grp_names)

  if len(btx) > 0:
    bottom_text(btx, pos=[0.05, 0.08])

  return ax1

def rtofs_reg2Dmaps():
  """
    Regional 2D maps
    1/12.5-degree HYCOM global grid and MOM6
  """
  REGN = {
    "NAtl1" : {
      "xl1" : 2330,
      "xl2" : 3340,
      "yl1" : 1570,
      "yl2" : 2090,
      "Reg" : 'NAtl',
      "ij1" : [2675, 1675],
      "ij2" : [2501, 1820],
      "ij3" : [2720, 1700],
      "Name": 'Trop/Subrtop N Atlantic'
    },
    "GINSea": {
      "xl1" : 2712,
      "xl2" : 3990,
      "yl1" : 2296,
      "yl2" : 3200,
      "Reg" : 'NPac',
      "ij1" : [2712, 2296],
      "ij2" : [3990, 3200],
      "ij3" : [3990, 2296],
      "Name": 'GIN Seas',
      "Ip"  : [2712, 2712, 3990, 3990],
      "Jp"  : [2296, 3200, 3200, 2296]
    },
    "Arctic": {
      "xl1" : 0,
      "xl2" : 4499,
      "yl1" : 2120,
      "yl2" : 3297,
      "Reg" : 'NPac',
      "ij1" : [1757, 2431],
      "ij2" : [1250, 2415],
      "ij3" : [1902, 2210],
      "Name": 'Arctic Ocean',
      "Ip"  : [352, 874, 1009, 1364, 1596, 1912, 2583, 2718, 3080, \
               3325, 3550, 3970, 4208],
      "Jp"  : [3286, 2655, 2591, 2606, 2654, 3265, 3281, 3084, 3056, 3013, \
               3008, 2644, 3287]
    },
    "Antrct": {
      "xl1" : 0,
      "xl2" : 4499,
      "yl1" : 0,
      "yl2" : 780,
      "Reg" : 'NPac',
      "ij1" : [1757, 2431],
      "ij2" : [1250, 2415],
      "ij3" : [1902, 2210],
      "Name": 'Southern Ocean',
      "Ip"  : [0, 0, 4499, 4499],
      "Jp"  : [0, 780, 780, 0]
    },
    "Agulhas": {
      "xl1" : 3400,
      "xl2" : 4000,
      "yl1" : 780,
      "yl2" : 1350,
      "Reg" : 'SAtl',
      "ij1" : [1757, 2431],
      "ij2" : [1250, 2415],
      "ij3" : [1902, 2210],
      "Name": 'Agulhas Region',
      "Ip"  : [3400, 3400, 4000, 4000],
      "Jp"  : [780, 1350, 1350, 780]
    },
    "NPac1" : {
      "xl1" : 661,
      "xl2" : 2028,
      "yl1" : 1800,
      "yl2" : 2612,
      "Reg" : 'NPac',
      "ij1" : [1757, 2431],
      "ij2" : [1250, 2415],
      "ij3" : [1902, 2210],
      "Name": 'Subpolar N Pacific',
      "Ip"  : [661,661,2028,2028],
      "Jp"  : [1793, 2355, 2355, 1793]
    },
    "SPac1" : {
      "xl1" : 29, 
      "xl2" : 1157,
      "yl1" : 1256,
      "yl2" : 1841,
      "Reg" : 'SPac',
      "ij1" : [648,  1442],
      "ij2" : [790,  1556],
      "ij3" : [964,  1332],
      "Name": 'Indonesia'
    },
    "GOM"     : {
      "xl1"   : 2300,
      "xl2"   : 2622,
      "yl1"   : 1690,
      "yl2"   : 1912,
      "Reg"   : 'NAtl',
      "Name"  : 'Gulf of Mexico',
      "Ip"    : [2425, 2333, 2326, 2545, 2567, 2573, 2546, 2518, 2482, 2450],
      "Jp"    : [1711, 1743, 1905, 1905, 1832, 1791, 1798, 1785, 1773, 1719]
    },
    "GulfStr" : {
      "xl1"   : 2512,
      "xl2"   : 3215,
      "yl1"   : 1793,
      "yl2"   : 2355,
      "Reg"   : 'NAtl',
      "Name"  : 'GulfStream',
      "Ip"    : [2512, 2512, 3215, 3215],
      "Jp"    : [1793, 2355, 2355, 1793]
    },
  }

  return REGN

def rtofs_sections():
  """
    Transect coordinates to plot vertical xsections
  """
  XSCT = {
    "GoM1"  : {
      "x1"  : 2340,
      "x2"  : 2560,
      "y0"  : 1798,
      "y1"  : 1700,
      "y2"  : 1915,
      "x0"  : 2503,
      "smin": 34.85,
      "smax": 36.75,
      "tmin": 3.85,
      "tmax": 27.85,
      "Reg" : 'NAtl',
      "Name": 'Gulf of Mexico'
    },
    "Portg" : {
      "x1"  : 3000,
      "x2"  : 3500,
      "y0"  : 2050,
      "y1"  : 1800,
      "y2"  : 2300,
      "x0"  : 3350,
      "smin": 34.92,
      "smax": 36.82,
      "tmin": 2.0,
      "tmax": 22.0,
      "Reg" : 'NAtl',
      "Name": 'East N Atlantic'
    },
    "NAtl1"  : {
      "x1"  : 2550,
      "x2"  : 3450,
      "y0"  : 1897,
      "y1"  : 1400,
      "y2"  : 2600,
      "x0"  : 3245,
      "Reg" : 'NAtl',
      "Name": 'Subtrop/Centr N Atlantic'
    },
    "GfStr" : {
      "x1"  : 2610,
      "x2"  : 2835,
      "y0"  : 1980,
      "y1"  : 1766,
      "y2"  : 2082,
      "x0"  : 2685,
      "Reg" : 'NAtl',
      "Name": 'West N Atlantic'
    },
    "GINS"  : {
      "x1"  : 3338,
      "x2"  : 3765,
      "y0"  : 2840,
      "y1"  : 2430,
      "y2"  : 3271,
      "x0"  : 3450,
      "smin": 33.50,
      "smax": 35.00,
      "tmin": -1.5,
      "tmax": 8.5,
      "Reg" : 'NAtl',
      "Name": 'GIN Seas'
    },
    "SAtl1" : {
      "x1"  : 3131,
      "x2"  : 3751,
      "y0"  : 1410,
      "y1"  : 920, 
      "y2"  : 1686,
      "x0"  : 3330,
      "smin": 34.37,
      "smax": 36.98,
      "tmin": 1.0,
      "tmax": 28.0,
      "Reg" : 'SAtl',
      "Name": 'Subtrop/Centr S Atlantic'
    },
    "SAtl2" : {
      "x1"  : 2706,
      "x2"  : 3970,
      "y0"  : 801, 
      "y1"  : 16,  
      "y2"  : 1220,
      "x0"  : 2050,
      "Reg" : 'SAtl',
      "Name": 'Southern S Atlantic'
    },
    "NPac1" : {
      "x1"  : 800,
      "x2"  : 1300,
      "y0"  : 2000,
      "y1"  : 1700,
      "y2"  : 2300,
      "x0"  : 1050,
      "Reg" : 'NPac',
      "Name": 'West N Pacific'
    },
    "NPac2" : {
      "x1"  : 887,
      "x2"  : 2022,
      "y0"  : 2221,
      "y1"  : 1948,
      "y2"  : 2617,
      "x0"  : 1390,
      "Reg" : 'NPac',
      "Name": 'Central N Pacific'
    },
    "NPac3" : {
      "x1"  : 809,
      "x2"  : 1667,
      "y0"  : 2995,
      "y1"  : 2649,
      "y2"  : 3286,
      "x0"  : 1325,
      "Reg" : 'NPac',
      "Name": 'Pacific Arctic Ocean'
    },
    "SPac1" : {
      "x1"  : 1510,
      "x2"  : 2580,
      "y0"  : 1000,
      "y1"  : 200,
      "y2"  : 1300,
      "x0"  : 1800,
      "Reg" : 'SPac',
      "Name": 'East S Pacific'
    },
    "SPac2" : {
      "x1"  : 1170,
      "x2"  : 2700,
      "y0"  : 575,
      "y1"  : 138,
      "y2"  : 1550,
      "x0"  : 1870,
      "Reg" : 'SPac',
      "Name": 'South/Centr S Pacific'
    },
    "IndO1" : {
      "x1"  : 4067,
      "x2"  : 796,
      "y0"  : 1423,
      "y1"  : 396, 
      "y2"  : 1500,
      "x0"  : 390, 
      "Reg" : 'IndO',
      "Name": 'Centr/East Indian Ocean'
    },
    "NPac4" : {
      "x1"  : 400,
      "x2"  : 700,
      "y0"  : 1601,
      "y1"  : 1500,
      "y2"  : 1700,
      "x0"  : 583,
      "Reg" : 'NPac',
      "Name": 'W. Pacific Polynesia'
    },
    "NPac5" : {
      "x1"  : 540,
      "x2"  : 620,
      "y0"  : 1608,
      "y1"  : 1595,
      "y2"  : 1640,
      "x0"  : 601,
      "Reg" : 'NPac',
      "Name": 'Sulu Sea'
    },
  }

  return XSCT

def colormap_salin(nclrs=200, clr_ramp=[1,1,1]):
  """
    Colormap for salinity
    low S value ramp to clr_ramp
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap
  import mod_colormaps as mclrs

  btm = cm.get_cmap('rainbow',nclrs)
  ixtop  = round(nclrs*0.1)-1  
  clrbtm = btm(range(nclrs))
  chbtm  = np.zeros((ixtop,4))
#
# Add ramp colors at the bottom of clrbar
#  if add_btm == True:
# Add white at the beginning:
  cxbtm  = clrbtm[0,:]
  
  chbtm[:,3] = cxbtm[3]

  for ik in range(3):
    cc0 = clr_ramp[ik]
    chbtm[:,ik]  = np.linspace(cxbtm[ik],cc0,ixtop)

  chbtm = np.flip(chbtm, axis=0)
  clrbtm = np.insert(clrbtm,0,chbtm, axis=0)

# Add extra colors at the top for better representation of 
# high-S range
  CLR = [[204,   0,   0],
         [153,   0,   0],
         [153,  76,   0],
         [204, 102,   0],
         [255, 229, 192]]
  CLR = np.array(CLR)/255.
  CLR[np.where(CLR > 1.)] =  1.
  CMP = mclrs.create_colormap(CLR, ixtop, cmp_obj=False)
  clr_high = CMP[0,:]

  nclrs  = clrbtm.shape[0]
  clrtop = clrbtm[-1,:]
  chtop  = np.zeros((ixtop,4))
  chtop[:,3] = cxbtm[3]
  for ik in range(3):
    cc0 = clr_high[ik]
    chtop[:,ik] = np.linspace(clrtop[ik],cc0,ixtop)

# COmbine high S colors at the end of colormap
  clrbtm = np.append(clrbtm, chtop, axis=0)
  clrbtm = np.append(clrbtm, CMP, axis=0)

  newclrs = clrbtm
  newcmp  = ListedColormap(newclrs)

  return newcmp
    
def colormap_temp(nclrs=200, clr_ramp=[1,1,1], add_btm=True):
  """
    Colormap for temp
    low S value ramp to clr_ramp
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  btm = cm.get_cmap('jet',nclrs)
  ixtop  = round(nclrs*0.1)-1
  clrbtm = btm(range(nclrs))
  chbtm  = np.zeros((ixtop,4))
  if add_btm == True:
# Add white at the beginning:
    cxbtm  = clrbtm[0,:]
  else:
# Add white at the top
    ixtop  = round(nclrs*0.1)-1
    ixbtm  = nclrs-ixtop-1
    cxbtm  = clrbtm[ixbtm,:]

  chbtm[:,3] = cxbtm[3]

  for ik in range(3):
    cc0 = clr_ramp[ik]
    chbtm[:,ik]  = np.linspace(cxbtm[ik],cc0,ixtop)

  if add_btm:
    chbtm = np.flip(chbtm, axis=0)
    clrbtm = np.insert(clrbtm,0,chbtm, axis=0)
  else:
    clrbtm[ixbtm+1:nclrs,:] = chbtm

  newclrs = clrbtm
  newcmp  = ListedColormap(newclrs)

  return newcmp

def colormap_ssh(cpos='Oranges',cneg='Blues_r',nclrs=100, clr_ramp=[1,1,1]):
  """
    Create colormaps for showing positive/negative ranges
    clr_ramp - color in the middle of the colormap
    nclrs - # of colors in positve or negative segments
  """
  from matplotlib import cm
  from matplotlib.colors import ListedColormap, LinearSegmentedColormap

  btm = cm.get_cmap(cpos,nclrs)
  ixpos  = round(nclrs*0.1)-1
  clrpos = btm(range(nclrs))
  chpos  = np.zeros((ixpos,4))
  chneg  = np.zeros((ixpos,4))

# Positive colors:
# Add white at the beginning:
  fixclr = clrpos[ixpos,:]
  chpos[:,3] = fixclr[3]

  clr_ramp=[1,1,1]
  for ik in range(3):
    cc0 = clr_ramp[ik]
    chpos[:,ik]  = np.linspace(fixclr[ik],cc0,ixpos)

  chpos = np.flip(chpos, axis=0)
  clrpos[:ixpos,:] = chpos

# Negative ssh range:
# Add white at the top
  clr0    = cm.get_cmap(cneg,nclrs)  
  clrneg  = clr0(range(nclrs))
  ixneg   = nclrs - ixpos - 1
  fixnegc = clrneg[ixneg,:]
  chneg[:,3] = fixnegc[3]

  for ik in range(3):
    cc0 = clr_ramp[ik]
    chneg[:,ik]  = np.linspace(fixnegc[ik],cc0,ixpos)

#  chneg = np.flip(chneg, axis=0)
  clrneg[ixneg+1:nclrs,:] = chneg

  newclrs = np.append(clrneg,clrpos, axis=0)
#  newclrs = clrbtm
  newcmp  = ListedColormap(newclrs)

  return newcmp
  

def plot_orthographic_proj(LON, LAT, A2d, cmpr=[], sttl='Orthographic Porjection',\
                           lon0=-10, lat0=40, res='l', rmin=[], rmax=[], \
                           btx=[], fgnmb=1):
  """
    Plot 2D maps in orthographic projections
    To change projection angle: set lon0, lat0
  """
  from mpl_toolkits.basemap import Basemap, cm
  import matplotlib.colors as colors
  import matplotlib.mlab as mlab
#  from matplotlib.colors import ListedColormap
  m = Basemap(projection='ortho', lon_0=lon0, lat_0=lat0, resolution=res)
#  m.drawcoastlines()
  ny = A2d.shape[0]; nx = A2d.shape[1]

# Check if colormap object is provided:
  try:
    Nclr = cmpr.colors
  except:
    cmpr = colormap_temp(clr_ramp=[0.9,0.8,1])
    cmpr.set_bad(color=[0.2, 0.2, 0.2])
#    cmpr.set_under(color=[1, 1, 1])  #
#  x, y = m(lons, lats) # compute map proj coordinates.
#  LONA = LON.copy()
#  LONA = np.append(LONA, np.array([LONA[:,0]]).transpose(), axis=1)
#  LONA = np.append(LONA, np.array([LONA[0,:]]), axis=0)
  xh, yh = m(LON,LAT)  # hycom coordinates on the projections

# Need to mask values that are not in the projection
  PMsk = ( (xh > 1e20) | (yh > 1e20) )
#  A2d[PMsk] = np.nan
  app = A2d.copy()
  app[PMsk] = np.nan

  xh[PMsk]  = 1.e30
  yh[PMsk]  = 1.e30

  if (not rmax) or (not rmin): 
    rmin, rmax = minmax_clrmap(app)

  fig1 = plt.figure(fgnmb,figsize=(9,9))
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])

#  im1 = m.pcolormesh(xh, yh, app)
#  im1.set_cmap(cmpr)
  im1 = m.pcolormesh(xh, yh, app, shading='nearest', cmap=cmpr,\
                     vmin=rmin, vmax=rmax)
# Plot land:
  PLnd = ( (~PMsk) & (np.isnan(A2d)) )
  xL = xh[PLnd]
  yL = yh[PLnd]
  Lclr = np.array([[0.2, 0.2, 0.2]])
  ax1.scatter(xL, yL, s=0.008, c=Lclr)
  

# alternatively but slower:
#  m.scatter(xh, yh, s=0.1, c=app, cmap=cmpr)
  m.drawparallels(np.arange(-90.,120.,10.))
  m.drawmeridians(np.arange(-180.,180.,10.))

  plt.title(sttl)

  ax2 = fig1.add_axes([ax1.get_position().x1+0.02,
               ax1.get_position().y0,0.02,
               ax1.get_position().height])
  clb = plt.colorbar(im1, cax=ax2, extend='both')
  ax2.set_yticklabels(ax2.get_yticks())
  ticklabs = clb.ax.get_yticklabels()
#  clb.ax.set_yticklabels(ticklabs,fontsize=10)
  clb.ax.set_yticklabels(["{:.2f}".format(i) for i in clb.get_ticks()], fontsize=10)
  clb.ax.tick_params(direction='in', length=5)

  if len(btx) > 0:
    bottom_text(btx, pos=[0.08, 0.02])

  plt.show()

  return
 

def minmax_clrmap(dmm, pmin=10, pmax=90, cpnt=0.01, fsym=False):
  """
  Find min/max limits for colormap 
  discarding pmin and 1-pmax min/max values
  cpnt - decimals to leave
  """
  dmm = dmm[~np.isnan(dmm)]
  a1  = np.percentile(dmm,pmin)
  a2  = np.percentile(dmm,pmax)
  cff = 1./cpnt
  rmin = cpnt*(int(a1*cff))
  rmax = cpnt*(int(a2*cff))

  if fsym and (rmin<0. and rmax>0.) :
    dmm = max([abs(rmin),abs(rmax)])
    rmin = -dmm
    rmax = dmm

  return rmin,rmax


def colormap_salin2(nclrs=200):
  """
    Colormap for salinity
  """
  import mod_colormaps as mclrs
  CLR = [[255, 255, 255],
         [255, 153, 204],
         [255, 102, 255],
         [255,   0, 255],
         [127,   0, 255],
         [178, 102, 255],
         [204, 153, 255],
         [153, 153, 255],
         [102, 102, 255],
         [0,   0,   255],
         [0,   102, 102],
         [51,  153, 255],
         [153, 204, 255],
         [153, 255, 204],
         [51,  255, 153],
         [0,   153,  76],
         [0,   204,   0],
         [102, 255, 102],
         [155, 255, 155],
         [255, 255, 153],
         [255, 255,   0],
         [204, 204,   0],
         [255, 128,   0],
         [255, 178, 102],
         [255, 204, 153],
         [255, 255, 204],
         [255, 153, 153],
         [255,  51,  51]]

  CLR = np.array(CLR)/255.
  CMP = mclrs.create_colormap(CLR, nclrs)

  return CMP



