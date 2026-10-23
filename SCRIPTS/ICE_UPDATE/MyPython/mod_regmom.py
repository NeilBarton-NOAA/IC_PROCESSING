"""
  Utility for regional MOM6-SIS2 
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import sys
#import pdb
import importlib
#import struct
from netCDF4 import Dataset as ncFile
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab
from matplotlib.patches import Polygon
from matplotlib.colors import ListedColormap
from mod_utils_fig import bottom_text
import mod_misc1 as mmisc1
import mod_bilinear as mblnr

def find_gridpnts_box(x0, y0, LON0, LAT, dhstep=0.5, \
                      ignore_north_lim=False, 
                      wrap_long=False,
                      use_close_indx=False):
  """
    Given pnt (x0,y0) find 4 grid points enclosing the pnt
    on a grid XX, YY - coordinates

     * -----------------  *
     |                    |
     |                    |
     |    *(x0,y0)        |
     |                    |
     |                    |
     |                    |
     * -----------------  *

    The algorithms works for both regular (straight orthogonal Mercator-type grids)
    and irregular (e.g. bipolar with slented grid lines) although in discontinuity regions
    where grid boxes are too narrow and slented may give errors 

    On a highly curvilinear portions of the domain (epsecially near the N. Pole)
    the algorithm may fail 
    because i-1, i+1, j-1, j+1 points may not be the closest points
    As a last option, allow to use 4 closest indices around the pnt
    use_close_indx=True <-- not recommended 


    Specify offset for i,j indices if a subset XX, YY is ued
    instead of whole domain X and Y (to speed up distance calculation)
    to get the global indices 

    dhstep - approximate max grid horizontal stepping (in units of XX/YY)
             in order to accelerate searching algorithm 
             making it too small may result in errors as
             grid points outside this range will be discarded 

    ignore_north_lim - True: 
           try to find points that are north of the northernmost
           grid points of the original grid
           This works for Polar stereogr. projection by grabbing points over the N. Pole 
           For other projections - this will likely not work

    wrap_long - True:
          when at the East/West bndry, search box vertices around the globe
          i.e. +/- 360 assuming grlobal grid, e.g. pnt at i=0 (0. deg) is connected to i=idm-1 (359.5 deg)

  """
  #import time
  #tt0 = time.perf_counter()
  import mod_bilinear as mblnr
  import mod_misc1 as mmisc1

# Need 2D arrays for LON, LAT
# if 1D array - Mercator grid is assumed
  ndim = len(LAT.shape)
  if ndim == 1:
    mm = len(LAT)
    nn = len(LON0)
    LON0 = np.tile(LON0,(mm,1))
    LAT  = np.tile(LAT,(nn,1)).transpose()

  mm,nn = LAT.shape

  # dhstep should be > max grid spacing in latitudes:
  max_dlat = max(
    np.nanmax(np.abs(np.diff(LAT, axis=0))),
    np.nanmax(np.abs(np.diff(LAT, axis=1)))
  )
  #print(max_dlat)

  assert dhstep > max_dlat, f"increase dhstep={dhstep} to be >= (lat grid dlt={max_dlat:.4f})"

  # Normalize x0
  x0 = (x0 + 360) % 360
  #tt1 = time.perf_counter()
  #print(f" dT1 = {tt1-tt0} sec")

  # Normalize LON0
  # Special treatment near discontinuities
  # Avoid +/- 180 and 0/360 discontinuites by shifting lon grid:
  if y0 < 90 - 2*max_dlat:
    LON = (LON0 + 360) % 360
  else:
    LON = mblnr.shift_longitudes(LON0, ref_lon=x0)  # expensive to run

  #tt2 = time.perf_counter()
  #print(f" dT2 = {tt2-tt1} sec")

  # Latitude range check - global min lat
  lat_min = np.nanmin(LAT)
  lat_max = np.nanmax(LAT)

  assert lat_min <= y0, f"check lat y0={y0:.4f} < min(LAT) {np.nanmin(LAT):.4f}"
  if not ignore_north_lim:
    assert lat_max >= y0, f"check lat y0={y0:.4f} > max(LAT) {np.nanmax(LAT):.4f}"

  # Check local min lat:
  # For curvilinear grids a point can be > global min(LAT) but still outside the domain
  # where the boundary "curves" when plotted in the lon/lat space
  dx = 2 * dhstep
  xl1 = x0 - dx
  xl2 = x0 + dx
  lon_mask = (LON > xl1)
  lon_mask &= (LON < xl2)
  if not lon_mask.any():
    raise AssertionError(f"Local min lat: No pnts found in long: {xl1:.4f}/{xl2:.4f}, pnt x0/y0: {x0:.4f}/{y0:.4f}")
   
  local_min = np.min(LAT[lon_mask]) 
  if y0 < local_min:
    print(f'pnt x0/y0: {x0:.3f}/{y0:.3f} outside: local min lat={np.min(LAT[JL,IL]):.4f} skipping ...')
    return [],[]

  # Subsample
  dy = dhstep

  # Away from Pole, limit Long. for quicker search
  if y0 < (90 - 2*max_dlat):
    # Wrapped longitude distance
    # Works for [0,360] !!!
    dlon = np.abs(LON - x0)
    dlon = np.minimum(dlon, 360.0 - dlon)  # this gives 1 - 359 = 2 not 358
    NPmask = (LAT > y0 - dy)
    NPmask &= (LAT < y0 + dy)
    NPmask &= (dlon < dx)
  else:
    # Keep all long. around N. Pole for wrapping over the N. Pole
    NPmask = (LAT > y0 - dy)
    NPmask &= (LAT < y0 + dy)

  JJ, II = np.where(NPmask)
  min_npnts = 8      # smaller N points - higher risk to not hav 4 pnts to include x0, y0
                     # near boundaries of tripolar grids
                     # too high N - slower code
  if len(JJ) < min_npnts:
    print(f"ERR: Subsampling the region around x0={x0}, y0={y0} failed")
    print(f"ERR: Only {len(JJ)} grid points are in ths subset")
    print(f"ERR: Try increasing dhstep={dhstep}")
    raise RuntimeError("Subsample region around x0,y0 failed: Not enough points")

  #tt3 = time.perf_counter()
  #print(f" dT3 = {tt3-tt2} sec")

  def find_closest_point(y0, x0, LON, LAT, JJ, II, Np=5):
    XX = LON[JJ, II]
    YY = LAT[JJ, II]
    DD = mmisc1.dist_sphcrd(y0, x0, YY, XX)
    kmin = np.argmin(DD)
    jmin = JJ[kmin]
    imin = II[kmin]
    xmin = LON[jmin, imin]
    ymin = LAT[jmin, imin]
    return jmin, imin, xmin, ymin

  def find_n_closest_points(y0, x0, LON, LAT, JJ, II, N=10):
    """
      On a curvilinear grid the closest grid vertex may belong to 
      cells that do not geometrically contain the target point.
      Try N closest point, if the 1st fails

      Better approach - find closest centroid of the grid boxes
    """
    XX = LON[JJ, II]
    YY = LAT[JJ, II]
    DD = mmisc1.dist_sphcrd(y0, x0, YY, XX)
    idx = np.argsort(DD)[:N]

    return JJ[idx], II[idx], DD[idx]  


  # Try finding N closest points and enclosing grid boxes for each of these
  # until find the right one
  # But better - find closest grid centroid not vertices
  JVX, IVX, _ = find_n_closest_points(y0, x0, LON, LAT, JJ, II, N=min_npnts)

  jv1, iv1 = JVX[0], IVX[0]  # keep this in case the find_box approach fails
  xv1 = LON[jv1,iv1]
  yv1 = LAT[jv1,iv1]

  #tt4 = time.perf_counter()
  #print(f" dT4 = {tt4-tt3} sec")

  INp = False
  if not INp:
    IV, JV, INp = find_box_include_comb(x0, y0, IVX, JVX, LON, LAT, eps_tol=1.e-2)

  # Not worked try i+/-1, j+/-1
  for jv, iv in zip(JVX, IVX):
    # skip boundaries
    if jv == 0:
      print(f'WARN: pnt x0/y0: {x0:.3f}/{y0:.3f} at the S boundary: i/j={iv1}/{jv1}, skipping ...')
      return [],[] 

    if jv == mm-1:
      if not ignore_north_lim:
        print(f'WARN: pnt x0/y0: {x0:.3f}/{y0:.3f} at N boundary: i/j={iv1}/{jv1}, skipping ...')
        return [],[] 
      else:
        # For grid with Merc. porjections, i.e. where N. Polar region
        # is split in halves and point on one side can be continued
        # to the other side over the N Boundary, take 4 closest pnts:
        ixx = IVX[:4]
        jxx = JVX[:4]    
        return ixx, jxx

    if not wrap_long and (iv == 0 or iv == nn-1):
      print(f'WARN: pnt x0/y0: {x0:.3f}/{y0:.3f} outside or near the E/W boundary: i/j={iv1}/{jv1}, skipping ...')
      return [],[] 

    if not INp:
      IV, JV, INp = find_box_include([x0,y0], [iv,jv], LON, LAT, eps_tol=1.e-8)

    if INp:
      #print("Found")
      #break
      ixx = np.array(IV).astype(int)
      jxx = np.array(JV).astype(int)
      return ixx, jxx

  # If nothing worked, try another approach
  # The following code is probably not needed:
  print(" --> Find_box methods failed, trying point search")

  # First guess for xv2: 
  # Find where the point lies wrt to the closest pnt:
  #                          * x(jv+1,iv1)
  #                          |
  #                          |
  #                          |
  #        *(x0,y0)          |
  #                          |
  #   *----------------------*-----------------* x(jv1,iv1+1)
  #  x(jv1,iv1-1)         x(jv1,iv1)
  #
  # Note that this approach will fail for curvilinear grid
  # when x0,y0 is very close to the closest point:
  #
  #                         x(jv1,iv1)
  #                         * -----------------* x(jv1,iv1+1)
  #                                                  
  #                         * x0,y0
  #   * (jv1,iv1-1)                                    
  
  # Select some reference point wrt to 1st closest pnt:
  xref, yref = xv1-0.01, yv1-0.01
  # Create projection centered at (xref,yref)
  transf_enu = mblnr.make_lonlat2xy_transformer(xref, yref)
  xv1c, yv1c = transf_enu.transform(xv1,yv1)
  x0c, y0c   = transf_enu.transform(x0,y0)
  #xv1c, yv1c = mblnr.lonlat2xy_enu(xv1, yv1, xref, yref) <-- identical transformatin but slower
  #x0c, y0c   = mblnr.lonlat2xy_enu(x0, y0, xref, yref)

  # Pnt orientation wrt I-axis:
  # Assumed that I is in the eastward direction
  # and locally xv1m < x0 < xv1p
  # Point at i-1:
  xv1m = LON[jv1,iv1-1]
  yv1m = LAT[jv1,iv1-1]
  xv1mc, yv1mc = transf_enu.transform(xv1m, yv1m)
  #xv1mc, yv1mc = mblnr.lonlat2xy_enu(xv1m,yv1m, xref, yref)

  # Point at i+1:
  xv1p = LON[jv1,iv1+1]
  yv1p = LAT[jv1,iv1+1]
  xv1pc, yv1pc = transf_enu.transform(xv1p,yv1p)
  #xv1pc, yv1pc = mblnr.lonlat2xy_enu(xv1p,yv1p, xref, yref)
  AI = np.array([xv1mc,yv1mc])   # pnt i-1
  BI = np.array([xv1pc,yv1pc])   # pnt i+1
  CP = np.array([xv1c,yv1c])     # closest pnt
  
  C = np.array([x0c,y0c])
  orntI = np.sign(mmisc1.orientation(AI,BI,C))
  if orntI == 0:
    orntI = 1     # pnt x0,y0 is exactly on I-axis
  # Check orientation of the closest pnt wrt to the x(i-1),x(i+1) line:
  orntI_CP = np.sign(mmisc1.orientation(AI,BI,CP))

  # Pnt orientation wrt J-axis:
  xv1j = LON[jv1-1,iv1]
  yv1j = LAT[jv1-1,iv1] 
  xv1jc, yv1jc = transf_enu.transform(xv1j,yv1j)
  #xv1jc, yv1jc = mblnr.lonlat2xy_enu(xv1j,yv1j, xref, yref)

  xv1l = LON[jv1+1,iv1]
  yv1l = LAT[jv1+1,iv1] 
  xv1lc, yv1lc = transf_enu.transform(xv1l,yv1l)
  #xv1lc, yv1lc = mblnr.lonlat2xy_enu(xv1l,yv1l, xref, yref)
  AJ = np.array([xv1jc,yv1jc])
  BJ = np.array([xv1lc,yv1lc])
  C = np.array([x0c,y0c])
  orntJ = np.sign(mmisc1.orientation(AJ,BJ,C))
  if orntJ == 0:
    orntJ = 1
  # Check orientation of the closest pnt wrt to the x(j-1),x(j+1) line:
  orntJ_CP = np.sign(mmisc1.orientation(AJ,BJ,CP))


  # Find pnt along I-axis in the same half-plane as x0,y0 wrt J-axis
  # Note the correct order of points when forming A, B
  C = np.array([xv1mc,yv1mc])
  orntIm1 = np.sign(mmisc1.orientation(AJ,BJ,C)) 
  C = np.array([xv1pc,yv1pc])
  orntIp1 = np.sign(mmisc1.orientation(AJ,BJ,C)) 

  if orntIm1 == 0 or orntIm1 == orntJ:
    iv2 = iv1-1
    jv2 = jv1
  elif orntIp1 == 0 or orntIp1 == orntJ:
    iv2 = iv1+1
    jv2 = jv1
  else:
    print(f"ERR: Could not find closest pnt 2 wrt to J-axis orientation, is x0,y0 outside domain?")
    raise Exception(f"x0={x0}, y0={y0}")

  xv2 = LON[jv2,iv2]
  yv2 = LAT[jv2,iv2]
  xv2c, yv2c = transf_enu.transform(xv2, yv2)
  #xv2c, yv2c = mblnr.lonlat2xy_enu(xv2, yv2, xref, yref)

  # In case pnt (x0,y0) lies exectly on the X-axis or Yaxis line between (xv1,yv1) and (xv2,yv2)
  # move it a bit to avoid possible errors in finding the other vertices:
  dxy = 1.e-9
  dist_pnt = mmisc1.distance_point2segment(x0,y0,xv1,yv1,xv2,yv1) 
  if dist_pnt < dxy:
    x0 = x0 + dxy
    y0 = y0 + dxy 
    x0c, y0c = transf_enu.transform(x0, y0)
    #x0c, y0c   = mblnr.lonlat2xy_enu(x0, y0, xref, yref)

# Define pnt orientation wrt to this line connecting 1st and 2nd vertex
#  to find which 
# direction to search for another vertex
# along J axis
# if at the domain bndry - ignore, nothing can do
  Dxv2 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv2c,yv2c],[x0c,y0c]))
  # If the pnt is on the X-line connecting [xv1c,yv1c],[xv2c,yv2c]
  #
  #            x0c,y0c
  #     *--------o--------*
  #
  if Dxv2 == 0.0:
    Dxv2 = 1.

  F_bndry = False
  if jv1+1 >= mm:
    Djp1 = - Dxv2
    F_bndry = True
  else:
    x1p1 = LON[jv1+1, iv1]
    y1p1 = LAT[jv1+1, iv1]
    x1p1c, y1p1c = transf_enu.transform(x1p1, y1p1)
    #x1p1c, y1p1c = mblnr.lonlat2xy_enu(x1p1, y1p1, xref, yref)
    Djp1 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv2c,yv2c],[x1p1c,y1p1c]))

# Singularity points:
# check if the grid points are same locations for varying j or i
# these points are over land - can be ignored
    D = mmisc1.dist_sphcrd(yv1, xv1, y1p1, x1p1)
    if D < 1.e-3:
      F_bndry = True
      Djp1 = - Dxv2
      
  if Djp1 != Dxv2:
# Try opposite direction from pnt 1 along J:
# Ignore the bndry
    if jv1-1 < 0:
      dirj = 1
      F_bndry = True
    else:
      x1m1 = LON[jv1-1, iv1]
      y1m1 = LAT[jv1-1, iv1]
      x1m1c, y1m1c = transf_enu.transform(x1m1, y1m1)
      #x1m1c, y1m1c = mblnr.lonlat2xy_enu(x1m1, y1m1, xref, yref)
      Djm1 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv2c,yv2c],[x1m1c,y1m1c]))

      if Djm1 != Dxv2 and not F_bndry:
        print(f"ERR: Cannot find direction to include pnt {x0} {y0}")
        raise Exception(f"Unexpected: orientations same Djp1={Djp1} Djm1={Djm1}")

      dirj = -1
  else:
    dirj = 1

  # Find a 3rd vertex along J-axis moving in direction = dirj 
  # such that the pnt (x0,y0) projected onto J falls within the 
  # segment (vrtx1 - vrtx3)
  iv3 = iv1
  jv3 = jv1 + dirj
  xv3 = LON[jv3,iv3]
  yv3 = LAT[jv3,iv3]
  xv3c, yv3c = transf_enu.transform(xv3, yv3)
  #xv3c, yv3c = mblnr.lonlat2xy_enu(xv3, yv3, xref, yref)
  Bv  = mmisc1.construct_vector([xv1c, yv1c],[xv3c, yv3c]) 
  Av  = mmisc1.construct_vector([xv1c, yv1c],[x0c, y0c])
  prA, cosT, tht = mmisc1.vector_projection(Av, Bv)
  lBv = np.sqrt(np.dot(Bv.transpose(), Bv))[0][0]
  prA = abs(prA)

  icc = 0
  while prA > lBv:
    # if projection is outside then extend the side to include the point in the box
    # may run ouside the domain!
    jv3 = jv3 + dirj
    if jv3 >= mm or jv3 < 0: 
      F_bndry = True
      jv3 = jv3 - dirj
      break

    xv3 = LON[jv3,iv3]
    yv3 = LAT[jv3,iv3]
    xv3c, yv3c = transf_enu.transform(xv3, yv3)
    #xv3c, yv3c = mblnr.lonlat2xy_enu(xv3, yv3, xref, yref)
    Bv  = mmisc1.construct_vector([xv1c, yv1c],[xv3c, yv3c]) 
    prA, cosT, tht = mmisc1.vector_projection(Av, Bv)
    lBv_new = np.sqrt(np.dot(Bv.transpose(), Bv))[0][0]
    if lBv_new < lBv:
      print(f"ERR: unexpected side xv3 not increasing, new: {lBv_new} old: {lBv}")
      raise Exception("unexpected edge length decreasing when choosing xv3")

    icc += 1
    if icc > 10:
      raise Exception(f"cannot find jv3 # steps={icc}")    

#ax1.plot(x0c,y0c,'r*')
#ax1.plot(xv1c,yv1c,'ro')
#ax1.plot(xv2c,yv2c,'go')
#ax1.plot(xv3c,yv3c,'bo')

# Check if first guess of xv2, yv2 was correct to include the point x0, y0:
  Dv2 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv3c,yv3c],[xv2c,yv2c]))
  Dp0 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv3c,yv3c],[x0c,y0c]))
  if Dp0 == 0: Dp0 = Dv2 # allow pnt to be on the line
  if Dv2 != Dp0:
# need to switch vertex x2 to opposite
# only one of those should be not 0:
    diri = np.sign(iv1-iv2)
    dirj = np.sign(jv1-jv2)

    icc = 0    
    while Dv2 != Dp0:
      iv2 = int(iv2 + diri)
      jv2 = int(jv2 + dirj)
      if iv2 < 0 or iv2 >= nn or jv2 < 0 or jv2 >= mm:
        F_bndry = True
        break
      xv2 = LON[jv2,iv2]
      yv2 = LAT[jv2,iv2]
      xv2c, yv2c = transf_enu.transform(xv2, yv2)
      #xv2c, yv2c = mblnr.lonlat2xy_enu(xv2, yv2, xref, yref)
      Dv2 = np.sign(mmisc1.orientation([xv1c,yv1c],[xv3c,yv3c],[xv2c,yv2c]))
  
      icc += 1 
      if icc > 10:
        if use_close_indx:
          print(f'Could not adjust vx2,, x0={x0:.4f}E, y0={y0:.4f}N, use 4 closest grid pnts')
          IV = [iv1-1, iv1,   iv1+1, iv1]
          JV = [jv1,   jv1-1, jv1,   jv1+1]
          ixx = np.array(IV).astype(int)
          jxx = np.array(JV).astype(int)
          return ixx, jxx
        else:
          raise Exception(f"adjusting vx2: too many steps {icc}")

  # Vertx 4: best guess use vertex 2 and vertex 3:
  iv4 = iv2
  jv4 = jv3

  # Construct a box enclosing the point x0, y0 and check: 
  IV = [iv1, iv2, iv4, iv3]
  JV = [jv1, jv2, jv4, jv3]
  IV0 = IV.copy()
  JV0 = JV.copy()

  # Convert box vertices into Cartesian coord wrt a reference pnt 
  XX  = LON[JV,IV]
  YY  = LAT[JV,IV]
  XXc, YYc = mmisc1.polygon_centroid(XX, YY) 
  #XXc, YYc = XX[0], YY[0]
  transf_enu = mblnr.make_lonlat2xy_transformer(XXc,YYc)
  XV, YV   = transf_enu.transform(XX,YY)
  x0c, y0c = transf_enu.transform(XXc,YYc) 
  #XV, YV  = mblnr.lonlat2xy_enu(XX, YY, XXc, YYc)
  #x0c, y0c = mblnr.lonlat2xy_enu(x0,y0, XXc, YYc)
  INp     = mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV)
  INp2    = mmisc1.inpolygon_1pnt(x0, y0, XX, YY)

  if F_bndry:
    if not INp:
      print(f"WARNING: x0={x0:6.2f} y0={y0:6.2f} on " +\
            f"bndry or out of the domain, cannot enclose")
    ixx = np.array(IV).astype(int)
    jxx = np.array(JV).astype(int)
    return ixx, jxx

# Most of cases should be handled by the above algorithm
# for very non-orthogonal I/J grid lines near singularities in 
# tripolar grids the point may still be outside the box
# Adjust sides to include the pnt if it is not in the box
  SIDE2V = {
     "side1": [0, 1],
     "side2": [1, 2],
     "side3": [2, 3],
     "side4": [3, 0]
  }

  V2SIDE = {
     "vertex0": [4, 1],
     "vertex1": [1, 2],
     "vertex2": [2, 3],
     "vertex3": [3, 4]
  }
     
# Plot box:
# ax1.plot(XV,YV,'.-')    
# ax1.plot(0,0,'r*')    # pnt should be inside XV,YV
# Check the grid points, if the grid is not rectilinear a regular boxes may not work
# find centroid for a simple polygon, and check x0,y0 positioning compared to the centroid
# should be on the same side for all sides of the box
# At least one side should give different orientations of centroid and the point x0, y0
  icc = 0
  while not INp:
    #xC = 0.25*(np.sum(XV))
    #yC = 0.25*(np.sum(YV))
    xC, yC = mmisc1.polygon_centroid(XV, YV)
    for iside in range(1,5):
      k1, k2 = SIDE2V[f"side{iside}"]  
      x1 = XV[k1]
      x2 = XV[k2]
      y1 = YV[k1]
      y2 = YV[k2]
      iv1 = IV[k1]
      iv2 = IV[k2]
      jv1 = JV[k1]
      jv2 = JV[k2]
      xg1 = XX[k1]
      yg1 = YY[k1]
      xg2 = XX[k2]
      yg2 = YY[k2]

      Dcentr = np.sign(mmisc1.orientation([x1,y1],[x2,y2],[xC,yC]))
      Dpnt0  = np.sign(mmisc1.orientation([x1,y1],[x2,y2],[x0c,y0c]))
#      print(f"side {iside} Dcentr={Dcentr} Dpnt0={Dpnt0}")
      if Dpnt0 == 0: Dpnt0 = Dcentr  # Allow point to be on the edge
      if Dpnt0 == Dcentr: continue   # this side is ok
# Find which vertex to move by checking which one is opposite to x0,y0 and centroid
      if Dcentr != Dpnt0:
        dx1xc = np.sign(xC-x1)
        dx1x0 = np.sign(x0c-x1)
        dy1yc = np.sign(yC-y1)
        dy1y0 = np.sign(y0c-y1)
        dx2xc = np.sign(xC-x2)
        dx2x0 = np.sign(x0c-x2)
        dy2yc = np.sign(yC-y2)
        dy2y0 = np.sign(y0c-y2)
        div1 = djv1 = div2 = djv2 = 0
        if dx1xc != dx1x0: div1 = dx1x0
        if dx2xc != dx2x0: div2 = dx2x0
        if dy1yc != dy1y0: djv1 = dy1y0
        if dy2yc != dy2y0: djv2 = dy2y0

        icc = 0
        while Dcentr != Dpnt0:
          iv1 = int(iv1+div1)
          iv2 = int(iv2+div2)
          jv1 = int(jv1+djv1)
          jv2 = int(jv2+djv2)

          IV[k1] = iv1
          IV[k2] = iv2
          JV[k1] = jv1
          JV[k2] = jv2
          XX  = LON[JV,IV]
          YY  = LAT[JV,IV]
          xC, yC = mmisc1.polygon_centroid(XX, YY)
          transf_enu = mblnr.make_lonlat2xy_transformer(xC, yC)
          XV, YV   = transf_enu.transform(XX,YY)
          x0c, y0c = transf_enu.transform(x0,y0)  
          #XV, YV  = mblnr.lonlat2xy_enu(XX, YY, xC, yC)
          #x0c, y0c = mblnr.lonlat2xy_enu(x0,y0, xC, yC)
          INp     = mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV)

          if INp: break

          x1 = XV[k1]
          x2 = XV[k2]
          y1 = YV[k1]
          y2 = YV[k2]
          Dcentr = np.sign(mmisc1.orientation([x1,y1],[x2,y2],[xC,yC]))
          Dpnt0  = np.sign(mmisc1.orientation([x1,y1],[x2,y2],[x0c,y0c]))

          icc += 1
          if icc > 10: 
            if use_close_indx:
              print(f'Could not find 4pnts, x0={x0:.4f}E, y0={y0:.4f}N, use 4 closest grid pnts')
              IV = [iv1-1, iv1,   iv1+1, iv1]
              JV = [jv1,   jv1-1, jv1,   jv1+1]
              INp = True
              break
            else:
              raise Exception(f"Couldnot adjust side to include point #iter={icc}")

    
  ixx = np.array(IV).astype(int)
  jxx = np.array(JV).astype(int)
  
  return ixx, jxx


def find_box_include_comb(x0, y0, IVX, JVX, LON, LAT, eps_tol=1e-3, eps0=1e-4):
  """
   Try combinations of N selected closest points in IVX, JVX to find
   a box that includes x0, y0
   making eps_tol too small - risk of false acceptance not quadrilateral boxes
                  too big  - rejecting semi-quadrilateral 
  """
  from itertools import combinations

  AL=list(combinations(range(len(IVX)), 4))
  for cp in AL:
    IV = IVX[list(cp)]
    JV = JVX[list(cp)]
    XX = LON[JV, IV]
    YY = LAT[JV, IV]

    # center longitude around target to avoid discont. over the fold
    XX = mblnr.shift_longitudes(XX, ref_lon=x0)

    # reject degenerate
    area = mmisc1.polygon_area(XX, YY)
    if abs(area) < 1e-12:
      continue

    # local projection 
    XXc, YYc = XX.mean(), YY.mean()
    transf = mblnr.make_lonlat2xy_transformer(XXc, YYc)
    XV, YV = transf.transform(XX, YY)
    x0c, y0c = transf.transform(x0, y0)
    XV, YV, IDX = mmisc1.reorder_polygon(XV, YV, indx=True)

    # reject self-intersecting quads
    if not mblnr.strictly_convex_quad(XV, YV, eps_tol=eps_tol):
      continue

    # inclusion test
    if mmisc1.point_on_edge(x0c, y0c, XV, YV, tol=eps0):
      return IV[IDX], JV[IDX], True

    if mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV, eps0=eps0):
      return IV[IDX], JV[IDX], True

  return [], [], False

def find_box_include(XY0, IJ1, LON, LAT, eps_tol=1.e-8):
  """
    Find a grid cell that encloses a pnt XY0
    given the first nearst vertex 

    Box #5 - larger box that inclues smaller 4 boxes
  """
  import mod_bilinear as mblnr

  x0, y0   = XY0 
  iv1, jv1 = IJ1
  BX = np.array([
      [[iv1,   jv1],
       [iv1,   jv1-1],
       [iv1-1, jv1-1],
       [iv1-1, jv1]],

      [[iv1,   jv1],
       [iv1+1, jv1],
       [iv1+1, jv1-1],
       [iv1,   jv1-1]],

      [[iv1,   jv1],
       [iv1,   jv1+1],
       [iv1+1, jv1+1],
       [iv1+1, jv1]],

      [[iv1,   jv1],
       [iv1-1, jv1],
       [iv1-1, jv1+1],
       [iv1,   jv1+1]],

      [[iv1-1, jv1-1],
       [iv1-1, jv1+1],
       [iv1+1, jv1+1],
       [iv1+1, jv1-1]],
  ])

  def inside_box(ibox):
    IV = BX[ibox,:,0]
    JV = BX[ibox,:,1]
    XX = LON[JV,IV]
    YY = LAT[JV,IV]

    # Check that this is not a degenerative polygon
    # all vertices on 1 line, etc - polygon area = 0
    area = mmisc1.polygon_area(XX, YY)
    if abs(area) < 1e-12:
      return False, None, None

    #XXc, YYc = mmisc1.polygon_centroid(XX,YY)
    XXc, YYc = XX[0], YY[0]
    transf_enu = mblnr.make_lonlat2xy_transformer(XXc, YYc)
    XV, YV   = transf_enu.transform(XX, YY)
    x0c, y0c = transf_enu.transform(x0, y0)
    #XV, YV   = mblnr.lonlat2xy_enu(XX, YY, XXc, YYc)
    #x0c, y0c = mblnr.lonlat2xy_enu(x0,y0, XXc, YYc)

    XV, YV  = mmisc1.reorder_polygon(XV,YV)
    # Check if the point is on one of the edges, count it as in the box 
    if mmisc1.point_on_edge(x0c, y0c, XV, YV, tol=eps_tol): 
      return True, IV, JV

    INp     = mmisc1.inpolygon_1pnt(x0c, y0c, XV, YV, eps0=eps_tol)

    return INp, IV, JV

  # Check larger box if one of the smaller boxes contain pnt:
  INp_big,IV_big, JV_big = inside_box(4)

  # Small grid boxes:
  for ibox in range(4):
    INp, IV, JV = inside_box(ibox)
    #print(f"Inside: {INp}")
    if INp:
      return IV, JV, True

  # If point is in the large box, 
  if INp_big:
    return IV_big, JV_big, True

  # Found nothing
  return [], [], False  

#ax1.cla()     
#ax1.plot(XV,YV,'.-')
#ax1.plot(XV[0],YV[0],'o')
#ax1.plot(x0c,y0c,'o')
#ax1.plot([XV[0],XV[-1]],[YV[0],YV[-1]],'-')
#ax1.cla()
#ax1.plot(XX,YY,'.-')
#ax1.plot(x0,y0,'o')
# mm,nn = LON.shape
#IN,JN = np.meshgrid(np.arange(nn),np.arange(mm))
#ax1.plot(IN,JN,'y.')
#ax1.plot(iv1,jv1,'ro')


def fill_npole(A2d, HLON, HLAT, HH, Rpole = 2.5, bad_val = None, npnts_max=10):
  """
    Fill North Pole hole if needed
    Note: the N. Pole "hole" region should be either NaN or some other "bad value" (e.g. 999)
    do distinguish it from valid values to be used for filling the gap
    Using 0's is not recommended unless 0 cannot be a valid value in the field   

    Rpole - radius of search domain around the North pole to locate NaNs and not nans
            for interpolation
    bad_val - provide missing values if other than NaN otherwise N. Pole may not be detected
    npnts_max - max number of the closest points for averaging
  """
  import mod_utils as mutil
  import mod_bilinear as mblnr

  mm,nn = HLAT.shape
  Acopy = A2d.copy()

  print(f'Filling North Pole hole R={Rpole:.2f}')
  # N. Pole region:
  lat_npole = 90. - Rpole
  if np.max(HLAT) < lat_npole:
    print(f"fill_npole: NPole beyond domain lat_npole={lat_npole:.2f}N, max HLAT={np.max(HLAT):.2f}")
    return A2d

  if bad_val is not None:
    NPmask = (HLAT >= lat_npole) & (A2d == bad_val)
    A2d[NPmask] = np.nan

  # NPole mask:
  # = -1 - outside the NPole region
  # =  1 - valid values
  # =  0 - land 
  # =  9 - N Pole hole 
  npole_dom  = (HLAT >= lat_npole)
  npole_land = npole_dom & (HH >= 0)
  npole_hole = npole_dom & (~npole_land) & np.isnan(A2d)
  npole_data = npole_dom & (~npole_land) & np.isfinite(A2d)

  NPmask = np.full((mm, nn), -1, dtype=np.int8)
  NPmask[npole_data] = 1
  NPmask[npole_land] = 0
  NPmask[npole_hole] = 9  

  JJ, II = np.where( NPmask == 9 )
  if JJ.size == 0:
    print("North Pole hole not found, nothing to fix, check NPole = NaN or set bad_val")
    return Acopy

  Xdon = HLON[NPmask == 1]
  Ydon = HLAT[NPmask == 1]
  Adata = A2d[NPmask == 1]

  for ipp in range(II.size):
    ii0 = II[ipp]
    jj0 = JJ[ipp]
    x0  = HLON[jj0,ii0]
    y0  = HLAT[jj0,ii0]

    DIST = mmisc1.dist_sphcrd(y0, x0, Ydon, Xdon)
    assert np.min(DIST) > 0., "fill_npole: unexpected 0 dist for donor points outside NPole hole"

    # Find N closest points based on DIST:
    INDX = np.argpartition(DIST, npnts_max)[:npnts_max]

    WT = 1./(DIST[INDX] + 1.e-6) # to avoid very small dist near N pole
    WT = WT/np.sum(WT)
    assert abs(1.-np.sum(WT)) < 1.e-8, "Check weights WT"

    trgt = np.sum(WT*Adata[INDX])
    assert (trgt >= np.min(Adata)) and (trgt <= np.max(Adata)),\
    f"Filling npole error: i={ii0} j={jj0} filled={trgt:.4f} min/max ={np.min(Adata):.4f}/{np.max(Adata):.4f}"

    Acopy[jj0,ii0] = trgt

  return Acopy

def smooth_edges_arctic(A2d, HLON, HLAT, HH, hlat0=65, Rsearch=0.25, fill_land=True, extrp='box'):
  """
    Smooth shapr gradient at the edges of the data field
    moothly damp nonzero values to zero within the Arctic region
    (HLAT >= hlat0) using distance-based weighting (extrp='inv_dist') or
    box-averaging (extrp='box')

    Zeros north of hlat0 are treated as ramp targets.

    fill_land = True: Also extend values over land to avoid gaps near the coast

    Rsearch - radius (degrees) where values are searched for extrapolation / ramping
  """

  Aex = A2d.copy()
  Aex[np.isnan(Aex)] = 0.
  if not fill_land:
    Aex[HH >= 0] = np.nan

  # Points to be filled:
  missing_data = (HLAT >= hlat0) & (Aex == 0)

  valid_src = (HLAT >= hlat0) & (Aex >= 0) # allow no-snow values as source pnts 
  hlat_src = HLAT[valid_src]
  hlon_src = HLON[valid_src]
  data_src = Aex[valid_src]
  #JS,IS = np.where( valid_src )

  JJ, II = np.where( missing_data )
  if JJ.size == 0:
    print("No missed data, nothing to fix")
    return A2d

  Rsearch_m = Rsearch*111e3  

  # Sort points by latitude south->north for smooth ramping
  # n2s = np.argsort(HLAT[JJ, II])

  print(f"Extrapolating/ramping 2D field in Arctic to {hlat0:.2f}")
  icc = 0
  for ipp in range(JJ.size):
    icc += 1
    if icc % 5000 == 0:
      prct = icc / float(JJ.size) * 100.
      print(f"  {prct:.2f}% done ...")

    ii0 = II[ipp]
    jj0 = JJ[ipp]
    x0  = HLON[jj0,ii0]
    y0  = HLAT[jj0,ii0]

    dlat = Rsearch
    dlon = Rsearch / max(np.cos(np.deg2rad(y0)), 1e-3)
    srch_box = (np.abs(hlat_src - y0) <= dlat) & \
               (np.abs(hlon_src - x0) <= dlon)
    if not np.any(srch_box):
      continue

    DD = mmisc1.dist_sphcrd(y0, x0, hlat_src[srch_box], hlon_src[srch_box])
    valid_data = (DD > 0) & (DD <= Rsearch_m)
    if not np.any(valid_data):
      continue

    DIST = DD[valid_data]
    Adata = data_src[srch_box][valid_data]

    if extrp == 'inv_dist':
      WT = 1./(DIST + 1.e-6) # to avoid very small dist
      WT = WT/np.sum(WT)
      trgt = np.sum(WT * Adata)

    elif extrp == 'box':
      # box average
      Adata = data_src[srch_box]
      trgt = np.nanmean(Adata)

    if trgt < np.min(Adata) - 1e-6 or trgt > np.max(Adata) + 1e-6:
        print(f"Warning: filled value {trgt:.4f} out of original data range at i={ii0} j={jj0}")

    Aex[jj0,ii0] = trgt

  return Aex

 
def extrapolate_to_lat_arctic(A2d, HLON, HLAT, HH, hlat0=65, Npnts=5, Rsearch=20., fill_land=True):
  """
    Extrapolate and smoothly damp nonzero values to zero within the Arctic region
    (HLAT >= hlat0) using distance-based weighting 
    Values gradually decrease towards the lat=hlat0 from valid data region > hlat0

    Rsearch - to speedup distance computation, find closest points within this R (degrees)

    fill_land = True: Also extend values over land to avoid gaps near the coast
    Npnts - # of the closest data points to use for averaging

  """

  HLON = (HLON + 360.) % 360

  Aex = A2d.copy()
  Aex[np.isnan(Aex)] = 0.
  if not fill_land:
    Aex[HH >= 0] = np.nan

  # Points to be filled:
  missing_data = (HLAT >= hlat0) & (Aex == 0)

  valid_src = (HLAT >= hlat0) & (Aex > 0) # do not allow no-snow values as source pnts 
  hlat_src = HLAT[valid_src]
  hlon_src = HLON[valid_src]
  data_src = Aex[valid_src]
  #JS,IS = np.where( valid_src )

  JJ, II = np.where( missing_data )
  if JJ.size == 0:
    print("No missed data, nothing to fix")
    return A2d

  print(f"Extrapolating/ramping 2D field in Arctic to {hlat0:.2f}")
  icc = 0
  for ipp in range(JJ.size):
    icc += 1
    if icc % 5000 == 0:
      prct = icc / float(JJ.size) * 100.
      print(f"  {prct:.2f}% done ...")

    ii0 = II[ipp]
    jj0 = JJ[ipp]
    x0  = HLON[jj0,ii0]
    y0  = HLAT[jj0,ii0]

    dlat = Rsearch
    dlon = Rsearch / max(np.cos(np.deg2rad(y0)), 1e-3)
    dlon_raw = np.abs(hlon_src - x0)
    dlon_wrap = np.minimum(dlon_raw, 360. - dlon_raw)  # wrap around the 0/360 discont
    srch_box = (np.abs(hlat_src - y0) <= dlat) & (dlon_wrap <= dlon)

    npsrch = np.count_nonzero(srch_box)
    if npsrch < Npnts:
      print(f"Not enough data pnts for Rsearch={Rsearch:.1f}, found {npsrch}, i={ii0} j={jj0}")
      if npsrch < 1:
        continue

    #DD = mmisc1.dist_sphcrd(y0, x0, hlat_src, hlon_src)
    DD = mmisc1.dist_sphcrd(y0, x0, hlat_src[srch_box], hlon_src[srch_box])

    # Choose first N smallest distances for averaging:
    Idist = np.argpartition(DD, Npnts-1)[:Npnts]
    DIST = DD[Idist] * 1e-3    # m --> km
    Adata = data_src[srch_box][Idist]
    dist_avrg = np.mean(DIST)
    data_avrg = np.mean(Adata) 

    # add closest point on the hlat0 latitude
    dist2lat0 = mmisc1.dist_sphcrd(y0, x0, hlat0, x0) * 1e-3
    data_lat0 = 0.

    dist_tg = np.array([dist_avrg, dist2lat0])
    data_tg = np.array([data_avrg, data_lat0])

    WT = 1./(dist_tg + 1.e-6) # to avoid very small dist
    WT = WT/np.sum(WT)
    trgt = np.sum(WT * data_tg)

    if trgt < np.min(data_tg) - 1e-6 or trgt > np.max(data_tg) + 1e-6:
        print(f"WARN: filled val{trgt:.4f} out of data range: {np.min(data_tg):.4f}/{np.max(data_tg):.4f} i={ii0} j={jj0}")

    Aex[jj0,ii0] = trgt

  return Aex

 
def fill_land(aa1,aa2,aa3,aa4,HH,A3d,JJ,II,Jocn,Iocn):
  """
    Fill all nans in 1D arrays - land points
    at least 1 array should be ocean pnt
  """
  la1  = len(np.where(np.isnan(aa1))[0])
  la2  = len(np.where(np.isnan(aa2))[0])
  la3  = len(np.where(np.isnan(aa3))[0])
  la4  = len(np.where(np.isnan(aa4))[0])
  kdmh = len(aa1)

  LL  = np.array([la1,la2,la3,la4])
  iL  = np.where(LL == kdmh)[0]
  inL = np.where(LL < kdmh)[0] 

# This should not happen - no land points
# assumed at least 1 array is nans
  if len(iL) == 0:
    return aa1, aa2, aa3, aa4

  AP  = np.column_stack((aa1,aa2,aa3,aa4))
# All land - find closes ocean point
  if len(iL) == 4:
#    raise Exception("All 4 1D arrays - land points, check Land masks")
    for ild in range(4):
      il1 = II[ild]
      jl1 = JJ[ild]    
      D   = np.square((Iocn-float(il1))**2 + (Jocn-float(jl1))**2)
      ixx = np.argmin(D)
      ioc = Iocn[ixx]
      joc = Jocn[ixx]
      AP[:, ild] = A3d[:,joc,ioc]

  APm = np.nanmean(AP, axis=1) 
  AP[:,iL] = APm[:,None]
 
  return AP[:,0].squeeze(), AP[:,1].squeeze(),\
         AP[:,2].squeeze(), AP[:,3].squeeze() 

def box_averaging(A2d, HH, box_size=3, land_fill=False, \
                  land_value=0., iS=None, iE=None, jS=None, jE=None, \
                  pole_wrap = False, LAT=None, LON=None):
  """
    equal-weight box averaging
    for smoothing 2D fields
    filtering is within A2d[jS:jE+1,iS:iE+1]
    at the domain boundaries - 1-side averaging if outside data are not available

    if land_fill == True:
    Land values are set to land_value 

  box_size - Full box size (e.g., 3 x 3)

   pole_wrap: True - special treatment of the grid (Lambert projection)
                    where the top of the grid is cut across the polar region
                    neighboring grid pnts are searched over the cut for
                    smooth averaging over this cut (wrapping values over the cut)

             False - averaging is performed to the jE or last row without 
                     searching for neighboring points across the polar cut
                     This is the right option for the not polar grid or 
                     where there is no wrapping issues (like polar porjections)
                    
  """
  from scipy.ndimage import uniform_filter

  # Make odd:
  if box_size % 2 == 0: 
    box_size += 1  

  jdm, idm = HH.shape
  if iS is None: iS = 0
  if iE is None: iE = idm - 1
  if jS is None: jS = 0
  if jE is None: jE = jdm - 1

  print(f"Box averging: box_size={box_size} subdomain: j/i: {jS}:{jE}/{iS}:{iE}, pole_wrap={pole_wrap}")

  # lon / lat required for polar wrapping:
  if pole_wrap and (LAT is None or LON is None):
    raise Exception(f"box_averaging: for pole_wrap LON and LAT required")

  # No polar wrapping for non-polar grid
  if pole_wrap and np.nanmax(LAT[jE, :]) < 89.0:
    print(f"box_averaging: grid does not reach N Pole {np.max(LAT[jE,:]):.2f},  pole_wrapping disabled\n")
    pole_wrap = False

  A = A2d.copy()
  LMsk = HH >= 0
  if land_fill:
    A[LMsk] = land_value
  else:
    A[LMsk] = np.nan

  # Ignore nans:
  valid = np.isfinite(A).astype(float)  # valid values
  A0 = np.nan_to_num(A, nan=0.0)  # replace nans to be used by uniform filter 
  dx = box_size // 2
  dy = box_size // 2

  if not pole_wrap:
    # Get avrg over boxes, ignoring nans, as these = 0, 1/N^2 * sum(A)
    # Get avrg valid (=1) and nans (=0) over boxes = 1/N^2 * sum([0,1,1,...]), i.e. N valid points/ N^2
    avrg_A0  = uniform_filter(A0, size=box_size, mode='nearest')
    avrg_pnts = uniform_filter(valid, size=box_size, mode='nearest')

    # Note N^2 is cancelled for both when avrg_A0/avrg_pnts
    #, i.e. this is simply sum(values)/N valid pnts
    AF = np.full_like(avrg_A0, np.nan)
    np.divide(avrg_A0, avrg_pnts, out=AF, where=avrg_pnts > 0)
    #AF[avrg_pnts == 0] = np.nan

  else:
    # Polar wrapping
    # Lambert-like grid is assumed
    # where the N Pole is split in 2 halves such that going along a latitude (e.g. 85 N)
    # one will "jump" from one half of the map to another
    # Strategy: create ghost cells extending the grid beyond the last row
    # and populate with the values from the last rows in reversed order

    # Find the shift point  line of the grid: 
    # where topmost latitude has minimum in a "M" shape 
    # second deriv > 0 and maximum for this case because
    # the polar cut creates a sharp change in curvature
    d2l = np.diff(LAT[-1,:], n=2) 
    ishift = np.argmax(d2l) + 1
    if abs(ishift - idm // 2) > 1: 
      print(f"WARN: Check shift line indx={ishift}, expected:  i={idm // 2} +/- 1") # should be close to idm/2

    #A_ghost = np.roll(A0[-dy:, ::-1], shift=idm//2, axis=1)
    #valid_ghost = np.roll(valid[-dy:, ::-1], shift=idm//2, axis=1)
    Acup = A0[-dy:, :]
    A_ghost = np.fliplr(np.flipud(Acup))
    valid_ghost = np.fliplr(np.flipud(valid[-dy:,:]))

    # Add ghost cells:
    A0_ext     = np.vstack([A0, A_ghost])
    valid_ext  = np.vstack([valid, valid_ghost])

    avrg_A0   = uniform_filter(A0_ext, size=box_size, mode='constant', cval=0.0)
    avrg_pnts = uniform_filter(valid_ext, size=box_size, mode='constant', cval=0.0)

    AF_ext = np.full_like(avrg_A0, np.nan)
    np.divide(avrg_A0, avrg_pnts, out=AF_ext, where=avrg_pnts > 0)

    # Extract the original domain:
    AF = AF_ext[:jdm, :]

  # Apply only in requested subdomain 
  Aout = A2d.copy()
  Aout[jS:jE+1, iS:iE+1] = AF[jS:jE+1, iS:iE+1]

  if not land_fill:
    Aout[LMsk] = np.nan

  return Aout

def check_bottom(AA):
  """
    Make sure there are no NaNs in the 1D vertical data
    AA is 1D array
    If HYCOM land mask does not match MOM's --> all nans in the 1D profile
    from HYCOM
  """
  if AA[0] == np.nan:
    print(f"1D profile: all values are nans, Land mask mismatch")
    return 

  izb = np.argwhere(AA == np.nan)
  if len(izb) == 0:
    return

  AA[izb] = AA[min(izb)-1]

  return AA
  
def derive_TSprof_WOA23(seas, YR, Xp, Yp, grd=0.25, conv2pot=True):
  """
    Extract T & S profiles from WOA23
    for specified locations Xp, Yp
    Note: T - in situ T's
    conv2pot: Convert in situ T to potential 
    WOA:
    # season: 1-12 monthly, 13-winter (Jan-Mar), 14-spring (Apr-Jun), ...
    https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/bin/woa23.pl
    woa23_[DECA]_[v][tp]_[gr].nc - NetCDF format
    where:
    [DECA] - decade
    [v] - variable
    [tp] - time period
    [ft] - field type
    [gr] - grid
    e.g: woa23_95A4_t14_04.nc
  """
  import mod_swstate as msw

  if grd==0.25:
    cgrd=4
  woa='woa23'

  if YR > 1990 and YR < 2005:
    deca = "95A4"
  elif YR > 2005 and YR < 2014:
    deca = "A5B4"
  elif YR > 2014 and YR < 2023:
    deca = "B5C2"

  urlT = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/' + \
         f'temperature/netcdf/{deca}/0.25/'
  urlS = 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/woa23/DATA/' + \
         f'salinity/netcdf/{deca}/0.25/'

  tfnm=f'{woa}_{deca}_t{seas:02d}_{cgrd:02d}.nc'
  sfnm=f'{woa}_{deca}_s{seas:02d}_{cgrd:02d}.nc'
#  print("Reading {1} from {0}".format(furl,varnm))
  
  tvar='t_an'
  svar='s_an'

  furl = os.path.join(urlT, tfnm)
  print(f"Extracting T prof from WOA23 decadal clim season={seas}")
  print(furl)
  nc   = ncFile(furl)
  T3d  = nc.variables[tvar][:].data.squeeze()
  T3d  = np.where(T3d>1.e10, np.nan, T3d)
  kdm, jdm, idm  = T3d.shape
 
  furl = os.path.join(urlS, sfnm)
  print(f"Extracting S prof from WOA23 decadal clim season={seas}")
  print(furl)
  nc   = ncFile(furl)
  S3d  = nc.variables[svar][:].data.squeeze()
  S3d  = np.where(S3d>1.e10, np.nan, S3d)

  ZZ   = nc.variables['depth'][:].data.squeeze()
  ZZ   = -abs(ZZ)
  lat  = nc.variables['lat'][:].data.squeeze()
  lon  = nc.variables['lon'][:].data.squeeze() 
  dlat = np.max(np.diff(lat))
  dlon = np.max(np.diff(lon))
 
  Npnts = len(Xp)
  print(f"N original locations={Npnts}, conv to potent temp={conv2pot}")

  if max(lon) <= 180. and min(lon) < 0.:
    Xp = np.where(Xp > 180., Xp-360., Xp)
 
  Iwoa = []
  Jwoa = []
  Tprf = np.zeros((kdm,1))
  Sprf = np.zeros((kdm,1))
  for kk in range(Npnts):
    x0 = Xp[kk]
    y0 = Yp[kk]

    dx = np.sqrt((lon-x0)**2)
    i0 = np.where(dx == min(dx))[0][0]
    dy = np.sqrt((lat-y0)**2)
    j0 = np.where(dy == min(dy))[0][0]
    chck = np.sqrt((x0-lon[i0])**2 + (y0-lat[j0])**2)
    if chck > dlat:
      print(f"ERR: WOA lon=lon[i0] lat=lat[j0]")
      raise Exception(f'ERR: no closest WOA grid pnt for x={x0} y={y0} kk={kk}')

    # check if this point already processed:
    if len(Iwoa) > 0:
      dI = int(np.min(np.sqrt((Iwoa-i0)**2 + (Jwoa-j0)**2)))
      if dI == 0: continue

    Iwoa.append(i0)
    Jwoa.append(j0)
    irc = len(Iwoa)
    tt  = T3d[:,j0,i0].squeeze()
    ss  = S3d[:,j0,i0].squeeze()
    if irc == 1:
      Tprf[:,0] = tt
      Sprf[:,0] = ss
    else:
      tt   = np.expand_dims(tt, axis=1)
      Tprf = np.append(Tprf, tt, axis=1)
      ss   = np.expand_dims(ss, axis=1)
      Sprf = np.append(Sprf, ss, axis=1)
    
  if conv2pot:
    pr_ref = 0.
    lat_obs = lat[Jwoa]
    for klr in range(kdm):
      t1d = Tprf[klr,:]
      s1d = Sprf[klr,:]
      inn = len(np.where(~np.isnan(t1d))[0])
      if inn == 0: break
      z0  = ZZ[klr]
      pr_db, pr_pa = msw.sw_press(z0, lat_obs)         
      Tpot = msw.sw_ptmp(s1d, t1d, pr_db, pr_ref)
      Tprf[klr,:] = Tpot

  return Tprf, Sprf, ZZ

def insitu2pot_3D(T3d, S3d, ZZ, HLAT, z_ref=0, uref='m'):
  """
    Convert in situ T to potential with pressure 
      reference: z_ref either in m (depth) 
      or dbar (pressure)
    T, S, Z - 3D arrays, lat0 - local latitude (2D array)
    ZZ - can be 1D or 3D array
  """
  import mod_swstate as msw

  kdm, jdm, idm = T3d.shape
  if len(ZZ.shape) == 1:
    Z3d = np.tile(ZZ, idm*jdm).reshape((idm,jdm,kdm))
    Z3d   = np.transpose(Z3d, (2, 1, 0))
  elif len(ZZ.shape) == 3:
    Z3d = ZZ
  else:
    raise Excpetion('ZZ array should be either 1D or 3D')

  Zref = np.zeros((jdm,idm)) + z_ref # ?? 
  if uref == 'm':
    if abs(z_ref) < 1.e-3:
      prref_db = np.zeros((jdm,idm))
      prref_pa = np.zeros((jdm,idm))
    else:
      prref_db, prref_pa = msw.sw_press(Zref, HLAT)
  else:
    prref_db = np.zeros((jdm,idm))
    prref_pa = np.zeros((jdm,idm))
 
  Tpot = np.zeros((kdm,jdm,idm))
  for klr in range(kdm):
    print(f'Converting T in situ --> T pot, layer {klr}')
    temp = T3d[klr,:].squeeze()
    sal  = S3d[klr,:].squeeze()
    z0   = Z3d[klr,:].squeeze()
    pr_db, pr_pa = msw.sw_press(z0, HLAT)
    tp = msw.sw_ptmp(sal, temp, pr_db, prref_db)
    Tpot[klr,:] = tp

  return Tpot

def insitu2pot_2D(T2d, S2d, zz0, HLAT, z_ref=0, uref='m'):
  """
    Convert in situ T to potential with pressure 
      reference: z_ref either in m (depth) 
      or dbar (pressure)
    T, S - 2D arrays, lat0 - local latitude 
    HLAT - 2D latitudes
    zz0 - in situ depth, m or dbar
  """
  import mod_swstate as msw

  jdm, idm = T2d.shape
  Zref = np.zeros((jdm,idm)) + z_ref
 
  if uref == 'm':
    if abs(z_ref) < 1.e-3:
      prref_db = np.zeros((jdm,idm))
      prref_pa = np.zeros((jdm,idm))
    else:
      prref_db, prref_pa = msw.sw_press(Zref, HLAT)
  else:
    prref_db = np.zeros((jdm,idm))
    prref_pa = np.zeros((jdm,idm))
 
  Tpot = np.zeros((jdm,idm))
  pr_db, pr_pa = msw.sw_press(zz0, HLAT)
  Tpot = msw.sw_ptmp(S2d, T2d, pr_db, prref_db)

  return Tpot

def insitu2pot_1D(T1d, S1d, Z1d, lat0, z_ref=0, uref='m', printT=False):
  """
    Convert in situ T to potential with pressure reference: z_ref either in m (depth) 
    or dbar (pressure)
    T, S, Z - 1D arrays, 1 profile, lat0 - local latitude
  """
  import mod_swstate as msw

  if uref == 'm':
    prref_db, prref_pa = msw.sw_press(z_ref, lat0)
  else:
    prref_db = z_ref

  Tpot = np.zeros((len(T1d)))
  for klr in range(len(T1d)):
    temp = T1d[klr]
    sal  = S1d[klr]
    z0   = Z1d[klr]
    if np.isnan(temp): continue
    pr_db, pr_pa = msw.sw_press(z0, lat0)
    Tpot[klr] = msw.sw_ptmp(sal, temp, pr_db, prref_db)

  if printT:
    mmisc.print_3col(Z1d,T1d,Tpot)

  return Tpot
 
def pot2insitu_1D(T1d, S1d, Z1d, lat0, z_ref=0, uref='m', printT=False):
  """
    Convert potential T to insitu with pressure reference for T potential: 
    z_ref either in m (depth) 
    or dbar (pressure)
    Z1d - depths (m) for computing in situ T
    T, S, Z - 1D arrays, 1 profile, lat0 - local latitude
  """
  import mod_swstate as msw

  if uref == 'm':
    pr_db, _ = msw.sw_press(z_ref, lat0)
  else:
    pr_db = z_ref

  Tsitu = np.zeros((len(T1d)))
  for klr in range(len(T1d)):
    temp = T1d[klr]
    sal  = S1d[klr]
    z0   = Z1d[klr]
    if np.isnan(temp): continue
    prref_db, _ = msw.sw_press(z0, lat0)
    Tsitu[klr] = msw.sw_ptmp(sal, temp, pr_db, prref_db)

  if printT:
    mmisc.print_3col(Z1d,T1d,Tsitu)

  return Tsitu










