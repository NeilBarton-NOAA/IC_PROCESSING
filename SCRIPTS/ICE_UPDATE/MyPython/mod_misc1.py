import numpy as np
import importlib

def inpolygon(xq, yq, xv, yv):
  """ 
  Function similar to matlab inpolygon
  based on a code from internet stackoverflow 
  returns "in" indicating if the query points specified by xq and yq 
  are inside or on the edge of the polygon area defined by xv and yv.
  """
  from matplotlib import path
  shape = xq.shape
  xq = xq.reshape(-1)
  yq = yq.reshape(-1)
  xv = xv.reshape(-1)
  yv = yv.reshape(-1)
  q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
  p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])

# Alternative:
# shape = xq.shape
# q = np.column_stack((xq.flatten(), yq.flatten()))
# p = path.Path(np.column_stack((xv.flatten(), yv.flatten())))

  return p.contains_points(q).reshape(shape)

def inpolygon_v2(X,Y,Xv,Yv):
  """
  Similar to inpolygon 
  X, Y - 2D numpy arrays of points
  Xv, Yv - coord of polygon
  returns: Mask of 0 and 1 that indicates where
           X,Y points are outside/inside
           IP, JP - indices of points inside the polygone
  """
  from matplotlib.path import Path

#  X, Y = np.meshgrid(np.arange(IDM), np.arange(JDM))
  IDM  = X.shape[1]
  JDM  = Y.shape[0]

  X, Y = X.flatten(), Y.flatten()
  pnts = np.vstack((Y,X)).T    # concatenate and transpose
  JI   = np.vstack((Yv,Xv)).T
  PP   = Path(JI)  # polygon
  grd  = PP.contains_points(pnts)
  MSK  = grd.reshape(JDM,IDM)
  MSK  = np.where(MSK,1,0)
  JP,IP = np.where(MSK==1)

  return MSK, IP, JP

def reorder_polygon(xv, yv, indx=False):
  """
    Order vertices of the polygon so that 
    they do not zig-zag: self-intersect sides (butterfly)
  """
  cx = np.mean(xv)
  cy = np.mean(yv)
  angles = np.arctan2(yv - cy, xv - cx)
  idx = np.argsort(angles)

  if indx:
    return xv[idx], yv[idx], idx
  else:
    return xv[idx], yv[idx]

def polygon_centroid(XV,YV):
  """
    Compute the centroid (Cx, Cy) of a 2D polygon.
    XV, YV: np array of (x, y) coords, ordered around the polygon.
  """
  AA = 0
  Cx = 0
  Cy = 0
  XV = np.asarray(XV).flatten()
  YV = np.asarray(YV).flatten()
  nn = len(XV)

  for ii in range(nn):
    # Compute the signed area
    # i.e. 1/2 of the sum of cross-products from each vertex
    # If the vertices are ordered counter-clockwise: A>0,
    # clockwise: A<0
    x0, y0 = XV[ii], YV[ii]
    x1, y1 = XV[(ii+1)%nn], YV[(ii+1)%nn]

    cross = x0*y1 - x1*y0
    AA += cross
    Cx += (x0 + x1) * cross
    Cy += (y0 + y1) * cross

  assert AA != 0, f"ERR: polygon_centroid: area 0 {AA}"
  AA *= 0.5
  Cx /= (6 * AA)
  Cy /= (6 * AA)

  return Cx, Cy

def point_on_segment(xq, yq, x1, y1, x2, y2, eps0=1e-8):
  """
    Check if point (xq, yq) lies on segment (x1, y1)-(x2, y2)
    Rule: point on the segment if:
    cross product = 0
    pnt is within the segment bound (dot product > 0 and < segm^2)
  """
  # Vectors
  dx = x2 - x1
  dy = y2 - y1
  dxp = xq - x1
  dyp = yq - y1

  # segm length:
  seglen = np.sqrt(dx*dx + dy*dy)
  #print(f"seglen={seglen}")
  seglen = np.max([seglen, 1e-10])

  # colinearity check (cross product = 0)
  cross = (dx*dyp - dy*dxp)/seglen  # normalized
  #print(f"cross={cross:.4e}")
  if abs(cross) > eps0:
    return False

  # check projection lies within segment bounds (dot >= 0 and <= |segment|^2)
  dot = dx*dxp + dy*dyp
  #print(f"dot = {dot:.4e}, seglen2={seglen**2:.4e}")
  if dot < 0:
    return False

  seglen2 = seglen**2
  if dot > seglen2 + eps0:
    return False

  return True

def point_on_edge(x0, y0, XV, YV, tol=1e-9):
  for k in range(len(XV)):
    x1, y1 = XV[k],     YV[k]
    x2, y2 = XV[(k+1)%len(XV)], YV[(k+1)%len(YV)]

    # compute distance from point to segment
    dist = distance_point2segment(x0, y0, x1, y1, x2, y2)
    #print(f"Distance pnt to edge: {dist:.4e}")
    if dist < tol:
      return True

  return False

def distance_point2segment(px, py, x1, y1, x2, y2):
  """
  distance_point_segment
  Return the minimum distance from point (px, py)
  to the segment (x1, y1) -> (x2, y2).
  """
  # segment vector
  vx = x2 - x1
  vy = y2 - y1

  # vector from start to point
  wx = px - x1
  wy = py - y1
  
  # project w onto v, compute parameter t
  # t < 0  -> closest is endpoint 1
  # t > 1  -> closest is endpoint 2
  # 0 <= t <= 1 -> closest is interior point
  Lsegm2 = vx*vx + vy*vy
  if Lsegm2 == 0:
    # the segment is a point
    dist_p2sgm = np.sqrt(wx**2 + wy**2)
    return dist_p2sgm
    
  t = (wx*vx + wy*vy) / Lsegm2
  t = max(0.0, min(1.0, t))   # clamp to the segment normalized [0,1]
    
  # closest point on the segment
  cx = x1 + t * vx 
  cy = y1 + t * vy
    
  # distance from point to closest point
  dist_p2sgm = np.sqrt((px-cx)**2 + (py - cy)**2)
  return dist_p2sgm

def inpolygon_1pnt(xq, yq, XV, YV, eps0=1e-8):
  """
    Point-in-polygon test: True if inside OR on edge
    standard winding number method
    eps0 - allowed error around 0
  """
  XV = np.asarray(XV).flatten()
  YV = np.asarray(YV).flatten()
  nn = len(XV)

  # (1) Edge test
  for ii in range(nn):
    x1, y1 = XV[ii], YV[ii]
    x2, y2 = XV[(ii+1)%nn], YV[(ii+1)%nn]
    if point_on_segment(xq, yq, x1, y1, x2, y2, eps0=eps0):
      return True

  # (2) Winding number test
  wnmb = 0

  for ii in range(nn):
    x1, y1 = XV[ii], YV[ii]
    x2, y2 = XV[(ii+1)%nn], YV[(ii+1)%nn]
    
    # Vector polygon edge, normalized:
    ax = x2 - x1
    ay = y2 - y1

    # Vector to the point q:
    qx = xq - x1
    qy = yq - y1

    # Cross produc vect_a x vect_q:
    AxQ = ax*qy - ay*qx

    # Count crossings
    if y1 <= yq < y2:
      # upward crossing
      #print(f"Upward: ii={ii}, AxQ={AxQ:.4e}")
      if AxQ > 0:
        wnmb += 1
    elif y2 <= yq < y1:
      # downward crossing
      #print(f"Downward: ii={ii}, AxQ={AxQ:.4e}")
      if AxQ < 0:
        wnmb -= 1

    #print(f"wnmb={wnmb}")

  #print(f"Winding number={wnmb}")
  return wnmb != 0

def polygon_area(X, Y):
  """
  Compute signed area of a 2D polygon using the shoelace formula.
  X, Y : arrays of coordinates ordered around the polygon.
  Returns signed area (positive for CCW, negative for CW).
  """
  X = np.asarray(X).flatten()
  Y = np.asarray(Y).flatten()
  n = len(X)

  area = 0.0
  for i in range(n):
    x0, y0 = X[i], Y[i]
    x1, y1 = X[(i+1) % n], Y[(i+1) % n]
    area += x0 * y1 - x1 * y0

  return 0.5 * area

def inpolygon_1pnt_v0(xq, yq, xv, yv):
  """ 
  Is 1 point(xq,yq) inside a polygon?
  Function similar to matlab inpolygon
  based on a code from internet stackoverflow
  returns True -  indicating if the query points specified by xq and yq 
  are inside or on the edge of the polygon area defined by xv and yv.
  Not very robust when the point is on the edge, perhaps need to increase radius ?
  """
  from matplotlib import path

  xv = np.asarray(xv).reshape(-1)
  yv = np.asarray(yv).reshape(-1)
  polygon = path.Path(np.column_stack((xv, yv)))

  q = np.array([[xq, yq]])

  # contains_point returns False for edge points; workaround: use radius=1e-10
  return polygon.contains_point((xq, yq), radius=1e-10)

def rotate_vector(uin,vin,thtd):
  """
  Rotate vector U(uin,vin)
  by angle thtd - in degrees
  """
  tht = thtd*np.pi/180.
  R = [np.cos(tht), -np.sin(tht), np.sin(tht), np.cos(tht)] # rotates c/clkws if tht<0
  R = np.array(R).reshape(2,2)
  UV = np.array([uin,vin]).reshape(2,1)
  UVr = R.dot(UV)

  ur = UVr[0].item()
  vr = UVr[1].item()

  """
  nf = 3
  arrowprops = dict(color='darkorange', linewidth=2)
  ax = compass(uin, vin, arrowprops, nf)

  nf = 4
  arrowprops = dict(color='blue', linewidth=2)
  ax = compass(ur, vr, arrowprops, nf)
  """

  return ur, vr


def date_yearday(YR,MM,DD):
  """
    Return day of the year
  """
  import time
  import datetime

  timeR = datetime.datetime(YR,1,1)
  timeN = datetime.datetime(YR,MM,DD)
  yrday = (timeN-timeR).days + 1

  return yrday
  

def datenum(ldate0,ldate_ref=[1,1,1]):
  """
  Given list [YY,MM,DD] - current date 
  compute days wrt to reference date - optional
  Similar to matlab datenum
  """
  import datetime

  ll = len(ldate0)
  YR = ldate0[0]
  MM = ldate0[1]
  DD = ldate0[2]
  if ll > 3:
    HH = ldate0[3]
    MN = ldate0[4]
  else:
    HH = 0
    MN = 0

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HHr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HHr = 0
    MNr = 0


  time0 = datetime.datetime(YR,MM,DD,HH,MN,0)
  timeR = datetime.datetime(YRr,MMr,DDr,HHr,MNr,0)

  dnmb = (time0-timeR).days+1+(HH-HHr)/24.+(MN-MNr)/1440.

  return dnmb


def datevec(dnmb,ldate_ref=[1,1,1]):
  """
  For datenum computed wrt to reference date - see datenum
  convert datenum back to [YR,MM,DD,HH,MN]
  """
  import time
  import datetime

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HHr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HHr = 0
    MNr = 0
  
  timeR = datetime.datetime(YRr,MMr,DDr,HHr,MNr,0)
  dfrct = dnmb-np.floor(dnmb)
  if abs(dfrct) < 1.e-6:
    HH = 0
    MN = 0
  else:
    HH = int(np.floor(dfrct*24.))
    MN = int(np.floor(dfrct*1440.-HH*60.))

  ndays = int(np.floor(dnmb))-1
  time0 = timeR+datetime.timedelta(days=ndays, seconds=(HH*3600 + MN*60))
  YR = time0.year
  MM = time0.month
  MD = time0.day
  HH = time0.hour
  MN = time0.minute
 
  dvec = [YR,MM,MD,HH,MN] 

  return dvec


def datevec1D(dnmb,ldate_ref=[1,1,1]):
  """
  For datenum computed wrt to reference date - see datenum
  convert datenum back to [YR,MM,DD,HH,MN]
  Input is 1D numpy array

  """
  import time
  import datetime

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HHr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HHr = 0
    MNr = 0
  
  timeR = datetime.datetime(YRr,MMr,DDr,HHr,MNr,0)
  dfrct = dnmb-np.floor(dnmb)
  HHi   = np.floor(dfrct*24.).astype(int)
  MNi   = np.floor(dfrct*1440.-HHi*60.).astype(int)

  ndays = (np.floor(dnmb)-1).astype(int)

  YR = []
  MM = []
  MD = []
  HH = []
  MN = []
  for it in range(np.shape(ndays)[0]):
    time0 = timeR+datetime.timedelta(days=ndays.item(it), \
                   seconds=(HHi.item(it)*3600 + MNi.item(it)*60))
    YR.append(time0.year)
    MM.append(time0.month)
    MD.append(time0.day)
    HH.append(time0.hour)
    MN.append(time0.minute)

  YR = np.array(YR)
  MM = np.array(MM)
  MD = np.array(MD)
  HH = np.array(HH)
  MN = np.array(MN) 
  dvec = [YR,MM,MD,HH,MN] 

  return dvec


def datestr(dnmb,ldate_ref=[1,1,1]):
  """
  For datenum computed wrt to reference date - see datenum
  convert datenum back to [YR,MM,DD,HH,MN]
  print the date
  """
  import time
  import datetime

  lr = len(ldate_ref)
  YRr = ldate_ref[0]
  MMr = ldate_ref[1]
  DDr = ldate_ref[2]
  if lr > 3:
    HHr = ldate_ref[3]
    MNr = ldate_ref[4]
  else:
    HHr = 0
    MNr = 0

  timeR = datetime.datetime(YRr,MMr,DDr,HHr,MNr,0)
  dfrct = dnmb-np.floor(dnmb)
  if abs(dfrct) < 1.e-6:
    HH = 0
    MN = 0
  else:
    HH = int(np.floor(dfrct*24.))
    MN = int(np.floor(dfrct*1440.-HH*60.))

  ndays = int(np.floor(dnmb))-1
  time0 = timeR+datetime.timedelta(days=ndays, seconds=(HH*3600 + MN*60))

  dstr = time0.strftime('%Y/%m/%d %H:%M')

  return dstr


def dist_sphcrd(xla1,xlo1,xla2,xlo2, Req=6371.0e3, Rpl=6357.e3):
  """
# this procedure calculates the great-circle distance ("spherical") between two
# geographical locations on an ellipse using Lambert formula
#
# lat-lon coordinates with its appropiate trigonometric
# signs. 
# INPUT: xla1, xlo1 - first point coordinates (latitude, longitude)
#        xla2, xlo2 - second point or 1D or 2D arrays of coordinates
#
# Coordinates can be 1 pnt, 1D or 2D arrays both - gives pnt by pnt distance or
# 1 pnt vs 1D or 2D - gives 1 pnt vs all other pnts distances
# 
# all input coordinates are in DEGREES: latitude from 90 (N) to -90,
# longitudes: from -180 to 180 or 0 to 360,
#
# keep longitidues in similar range, -180/180 or 0/360
#
# in the latter case, distances from Pnt 1 (LAT1,LON1) to all pnts (LAT2,LON2)
# are calculated 
# OUTPUT - distance (in m)
# R of the earth is taken 6371.0 km
#
  """
# print("xla1=",xla1)
# breakpoint()
  xla1 = np.float64(xla1)
  xlo1 = np.float64(xlo1)
  xla2 = np.float64(xla2)
  xlo2 = np.float64(xlo2)

  flt1 = False
  flt2 = False
  if isinstance(xlo1, float):
    flt1 = True
  if isinstance(xlo2, float):
    flt2 = True
#  else:
#    raise Exception("Either 1st or 2nd lon/lat have to be a single value")

  if np.absolute(xla1).max() > 90.0:
    print(f"ERR: dist_sphcrd Lat1 > 90 {xla1}.max()")
    dist = float("nan")
    return dist
  if np.absolute(xla2).max() > 90.0:
    print(f"ERR: dist_sphcrd Lat2 > 90  {xla2}.max()")
    dist = float("nan")
    return dist

# Convert into -180/180 - works best
#  if flt1: 
#    if xlo1 > 180.:
#      xlo1 = xlo1-360.
#  else:
#    if np.max(xlo1) > 180.:
#      xlo1 = np.where(xlo1>180., xlo1-360., xlo1)
#
#  if flt2: 
#    if xlo2 > 180.:
#      xlo2 = xlo2-360.
#  else:
#    if np.max(xlo2) > 180.:
#      xlo2 = np.where(xlo2>180., xlo2-360., xlo2)

  cf = np.pi/180.
  phi1 = xla1*cf
  phi2 = xla2*cf
  lmb1 = xlo1*cf
  lmb2 = xlo2*cf
  dphi = abs(phi2-phi1)
  dlmb = abs(lmb2-lmb1)

  I0s  = []
  J0s  = []
  epsD = 1.e-12 
# if type(dphi) in (int,float):   # scalar input
  if np.min(dphi) < epsD or np.min(dlmb) < epsD:
    if flt1 and flt2:
      if dphi == 0.0 and dlmb == 0.0:
        dist_lmbrt = 0.0
        return dist_lmbrt
    else:           # N dim array
      if dphi.ndim == 1:
        J0s = np.where((dphi<epsD) & (dlmb<epsD))
        if len(J0s) > 0:
          dphi[J0s] = epsD
          dlmb[J0s] = epsD
      elif dphi.ndim == 2:
        J0s, I0s = np.where((dphi<epsD) & (dlmb<epsD))
        if len(I0s)>0:
          dphi[J0s,I0s] = epsD
          dlmb[J0s,I0s] = epsD
#
  Eflat = (Req-Rpl)/Req  # flatenning of the Earth
# Haversine formula to calculate central angle:
  aa1 = (np.sin(dphi/2.))**2
  aa2 = np.cos(phi1)*np.cos(phi2)*(np.sin(dlmb/2.))**2
  dsgm_hv = 2.*np.arcsin(np.sqrt(aa1+aa2))  # haversine form for central angle
#
# Reduced latitides - due to flattening:
  beta1 = np.arctan((1.-Eflat)*np.tan(phi1))
  beta2 = np.arctan((1.-Eflat)*np.tan(phi2))
  PP = 0.5*(beta1+beta2)
  QQ = 0.5*(beta2-beta1)
  X = (dsgm_hv-np.sin(dsgm_hv))*( (np.sin(PP))**2*(np.cos(QQ))**2 )/( (np.cos(dsgm_hv/2.))**2 )
  Y = (dsgm_hv+np.sin(dsgm_hv))*( (np.cos(PP))**2*(np.sin(QQ))**2 )/( (np.sin(dsgm_hv/2.))**2 )
# if np.sin(dsgm_hv/2.) == 0.0:
#  breakpoint()

  dist_lmbrt = Req*(dsgm_hv-Eflat/2.*(X+Y))

  if np.min(dist_lmbrt)<0.0:
    print('WARNING: spheric distance <0: ',np.min(dist_lmbrt))

#  print(f"Min dist = {np.min(dist_lmbrt)}")

# breakpoint()

  return dist_lmbrt


def box_fltr(AA, i1=-1, i2=-1, j1=-1, j2=-1, dist_wgt='linear',nbx=9):
  """
    Box filter
    sqrt of box-size should be odd = 3, 5, 7, 9, ... 
    so that the central node would have same size in all directions!

    i1,i2, j1,j2 - start/end of the subdomain or whole domain if i1 <0, ... 
    dist_wgt = 'equal' weights distributed equally 1/nbx
             = 'linear' weights linearly decrease from the central pnt
             = 'invdst' inverse distance decrease 
  """
  jdm = AA.shape[0]
  idm = AA.shape[1]
  if i1 < 0 or i2 < 0 or j1 < 0 or j2 < 0:
    i1 = 0
    i2 = idm-1
    j1 = 0
    j2 = jdm-1

# define weights:
  ibx = int(np.sqrt(nbx))
  if not ibx*ibx == nbx:
    raise Exception('SQRT of Filter box size should give integer nbx={0}'.\
        format(nbx))

# Check size of the box side should be odd:
  if int((ibx-1)/2)*2 == ibx-1:
    print('Box-smoothing {0} weights, nbx={1} ...'.format(dist_wgt,nbx))
  else:
    ibx = ibx+1
    print('Boxsize {0} has even dx/dy, adjusted to nbx={1}, start filtering ...'.\
         format(nbx, ibx*ibx))
    nbx = ibx*ibx

  jbx = ibx
  ibx0 = int((ibx-1)/2)    # -1 for python index: central node index
  WGT = np.zeros((ibx,ibx))
  if dist_wgt == 'equal':
    WGT = WGT + 1./float(nbx)
  elif dist_wgt == 'linear':
    R0 = np.sqrt(2.*(ibx0+0.5)**2)
    for iib in range(ibx):
      for jjb in range(jbx):
        rr = R0-np.sqrt(float(iib-ibx0)**2 + float(jjb-ibx0)**2)
        WGT[jjb,iib] = rr
    
  elif dist_wgt == 'invdst':
    for iib in range(ibx):
      for jjb in range(jbx):
        rr = np.sqrt(float(iib-ibx0)**2 + float(jjb-ibx0)**2)
        if rr < 1.e-12:
          r05 = np.sqrt(0.25**2)
          WGT[jjb,iib] = 1/r05
        else:
          WGT[jjb,iib] = 1./rr

# Normalize so that sum(WGT) = 1
  smm = np.sum(WGT)
  WGT = WGT/smm

  Aflt = AA.copy()
  dibx = int((ibx-1)/2)
  djbx = int((jbx-1)/2)
  for ii0 in range(i1,i2):
    for jj0 in range(j1,j2):
      ii1 = ii0-dibx
      ii2 = ii0+dibx
      jj1 = jj0-djbx
      jj2 = jj0+djbx

      if ii1 < 0 or ii2 == idm or jj1 < 0 or jj2 == jdm:
        continue

      Asub = AA[jj1:jj2+1, ii1:ii2+1]
      aa0 = np.sum(Asub*WGT)

      Aflt[jj0,ii0] = aa0


  return Aflt

def print_1col(A,wd=8,prc=2):
  """
    Print out 1 colume of data
  """
  A = np.asarray(A)

  if A.ndim != 1:
    raise ValueError("Input must be 1-dimensional")

  for i, value in enumerate(A, start=1):
    print(f"{i}: {value:{width}.{precision}f}")

  return

def print_2col(A1, A2, wd=8, prc=2, kend=None):
  """
    Print out 2 columns of data, same size or not

    If arrays have different lengths, stop at the shorter one.
    If kend is specified, stop at kend.
  """
  A1 = np.asarray(A1)
  A2 = np.asarray(A2)

  ndim1 = A1.shape[0]
  ndim2 = A2.shape[0]
  if kend is not None:
    kend = min([ndim1, ndim2])

  for k, (val1, val2) in enumerate(zip(A1[:kend], A2[:kend]), start=1):
    print(f"{k}: {val1:{width}.{precision}f}  {val2:{width}.{precision}f}")

  return

def print_cols(*arrays, wd=15, prc=8, kend=None, str=None):
  """
  Print N columns of numeric data.
  Provide N 1D arrays or lists:
  a = [1.1, 2.2, 3.3]
  b = [10.1, 20.2, 30.3]
  c = [100.1, 200.2, 300.3]
  ...

  Stops at the shortest array length unless kend is specified

  str - optional title, e.g. str='snon i=100 j=105'
  """
  arrays = [np.asarray(a) for a in arrays]

  n = min(len(a) for a in arrays)

  if kend is not None:
    n = min(n, kend)

  if str is not None:
    print(f"{str}")

  for k, vals in enumerate(zip(*(a[:n] for a in arrays)), start=1):
    s = "  ".join(f"{v:{wd}.{prc}f}" for v in vals)
    print(f"{k}: {s}")

  return


def polysegm_indx(IV, JV):
  """
    Find index points connecting verticies of 
    the polysegment curve

    * 1 IV[0], JV[0]

         * 2 IV[1], JV[1]

     * 3  IV[2], JV[2]

  """
  nseg = len(IV)-1
  II   = []
  JJ   = []

  for isgm in range(nseg):
    I1 = IV[isgm]
    J1 = JV[isgm]
    I2 = IV[isgm+1]
    J2 = JV[isgm+1]

    IS, JS = xsect_indx([I1,I2],[J1,J2])
    if isgm > 0:
      IS = IS[1:]
      JS = JS[1:]

    if isgm == 0:
      II = IS.copy()
      JJ = JS.copy()
    else:
      II = np.append(II,IS)
      JJ = np.append(JJ,JS)

#  II = np.array(II, dtype='int')
#  JJ = np.array(JJ, dtype='int')
  return II, JJ


def xsect_indx(IIv, JJv):
  """
    Simple segment: 2 end points only
    IIv=[Is,Ie], JJv=[Js,Je]
    Find index points connecting 
    2 vertices (is,js) and (ie,je)
    Finding the shortest distance
      
                  * (Ie, Je)
                  |
            *--*--*
            |
      *--*--*
      |
      *   (Is, Js)

  """
  Is, Ie = IIv
  Js, Je = JJv
  II = np.array([Is])
  JJ = np.array([Js])
  cntr = 0
  x0   = Is
  y0   = Js
  ss0  = np.sqrt((Is-Ie)**2 + (Js-Je)**2)
  ss   = 1.5*ss0

  while ss > 0.5:
    di   = (Ie-x0)/ss
    dj   = (Je-y0)/ss
    Xnrm = np.array([-dj,di])    # normal to the X-section
    ss   = np.sqrt((x0-Ie)**2+(y0-Je)**2)
    if ss < 0.5:
      break
    x0   = x0+di
    y0   = y0+dj
    if np.round(x0)==II[-1] and np.round(y0)==JJ[-1]:
      continue
   
# Make sure that di and dj < 1 to avoid skipping of grid point 
# at line zgzaging
    xnew = np.round(x0)
    ynew = np.round(y0)
    dx   = xnew - II[-1]
    dy   = ynew - JJ[-1]
    nl  = np.sqrt(dx**2+dy**2)
    if nl > 1.:
# Add missing point closest to the main section
# dx and dy should be close to 1
      if abs(np.round(dx)) > 1.0 or abs(np.round(dy)) > 1.0:
        raise Exception(f"Too large dx={dx} or dy={dy} > 1")

# Moving along X
      xpi = II[-1] + dx
      ypi = JJ[-1]
# Moving along Y:
      xpj = II[-1]
      ypj = JJ[-1] + dy
#
# Find close pnt to the section line choosing the smallest
# area of triangles btw the pnt and the line 
      TrI = np.array([[xpi, ypi, 1],[x0, y0, 1],[x0-di, y0-dj, 1]])
      TrJ = np.array([[xpj, ypj, 1],[x0, y0, 1],[x0-di, y0-dj, 1]])
      ArI = abs(0.5*np.linalg.det(TrI))
      ArJ = abs(0.5*np.linalg.det(TrJ))
      if ArI < ArJ:
        II = np.append(II, np.round(xpi))
        JJ = np.append(JJ, np.round(ypi))
      else: 
        II = np.append(II, np.round(xpj))
        JJ = np.append(JJ, np.round(ypj))

    cntr  += 1
    II   = np.append(II,np.round(x0))
    JJ   = np.append(JJ,np.round(y0))

    if cntr > 1e6:
      print(f"Couldnt find indices? Too many points: {cntr}")
      raise Exception("Finding index points: Endless loop ")

    if ss>2.*ss0:
      raise Exception("Finding index points: Moving away from the segment")
      
  dE = np.sqrt((II[-1]-Ie)**2 + (JJ[-1]-Je)**2)
  if dE > 1.e-6:
    II = np.append(II, Ie)
    JJ = np.append(JJ, Je)

  II = II.astype(int)
  JJ = JJ.astype(int)  
  return II, JJ

class SEGMENTS():
  def __init__(self, Isgm, Jsgm, Lsgm1, Lsgm2, vnrm1, vnrm2, nLeg, \
               curve_ornt, LegNorm):
    self.curve_ornt = curve_ornt
    self.I_indx     = Isgm
    self.J_indx     = Jsgm
    self.half_Lsgm1 = Lsgm1
    self.half_Lsgm2 = Lsgm2
    self.half_norm1 = vnrm1
    self.half_norm2 = vnrm2
    self.Leg_number = nLeg
    self.LegNorm   = LegNorm
#    self.Tht2Sct1   = Tht2Sct1
#    self.Tht2Sct2   = Tht2Sct2

def define_segments(IJ, DX, DY, curve_ornt='positive', check_pole=True):
  """
    Given poly-segment section connected at vertices IJ[iS:iE, jS:jE]
    find grid indices connecting the vertices

    For more accurate volume/ heat/ FW transport calculation:
    compute fluxes via halves of the grid cell allowing 
    "staircase" structure of the sections approximating
    slanted section that does not coincide with the model
    grid axes

    Each segment consists of 2 halvs wrt to the center pnt
    For each half norm is defined as a negatively oriented curve

           half2 of isgm
         *---------------------* (isgm+1)
         |     | nrm2
 half1   |-->  V    :
 of isgm | nrm1     :
         |----------:-------------------
         |          :
         |          :
isgm-1   *          :
         |          :
         |          :
         |----------:-------------------



    find segment lengths for each half
    
    define positive norm - depending on the curve or section
                           orientation following math definition:
     "A positively oriented curve is a planar simple closed curve 
     such that when traveling on it one always has the 
     curve interior to the left"

     Note right or left depends on the segment indices and 
     negative/positive is defined using given index order

    To preserve volume flux: U should be multiplied by the norm

    All small segments are combined into 1 array

    For plotting U sections, U-vector should be projected onto
    actual norm vector of the section line(s) 
    This norm is LegNorm (is the same for 1 leg) 
    
    check_pole - default True, meaning global grid has discontinuity
            at the N Pole, i.e. a line across the Arctic Ocean (Alaska - Fram)
            has 2 discontinuous segments on the grid

  """
  jdm      = DX.shape[0]
  idm      = DX.shape[1]
  nlegs    = IJ.shape[0]-1
  nLeg     = np.array([])
  vnrm1    = np.array([])  # norm vector half-segment 1
  vnrm2    = np.array([])
  Lsgm1    = np.array([])  # length of half-segment 1
  Lsgm2    = np.array([])
  Isgm     = np.array([])  # indices of segments along section 
  Jsgm     = np.array([])
#  Tht2Sct1 = np.array([]) # angle from half-segm to main section 
#  Tht2Sct2 = np.array([]) # sign - math. convention
  LegNorm  = np.array([]) # norm vector for the main sections
                          # normal orientation wrt to curve_ornt
  for ileg in range(nlegs):
    Is  = IJ[ileg,0]
    Js  = IJ[ileg,1]
    Ie  = IJ[ileg+1,0]
    Je  = IJ[ileg+1,1]
    IIv = [Is,Ie]
    JJv = [Js,Je]
#
# Check if this segment is over the N. Pole:
    if int(Js) == (jdm-1) and int(Je) == (jdm-1) and (check_pole):
      print(f'Segment is over the N. Pole, skipping leg {ileg}')
      continue

#    II, JJ = mmisc.xsect_indx(IIv,JJv)
    II, JJ = xsect_indx(IIv,JJv)  # index of grid points, not half-segm!

# Find norm vector for the main section
# Leg vector defining the main direction:
    LegV = np.array([IIv[1]-IIv[0],JJv[1]-JJv[0],0])
#    LegV = np.expand_dims(LegV, axis=0)
    LegV_mgn = np.sqrt(np.dot(LegV.transpose(),LegV))
    LegVnorm = find_norm([IIv[0],JJv[0]], [IIv[1],JJv[1]], curve_ornt)
    LegVnorm = np.expand_dims(LegVnorm, axis=0) 

    Isgm    = np.append(Isgm, II)
    Jsgm    = np.append(Jsgm, JJ)
#    print(f'(1) Isgm={Isgm[:10]}')

    nsgm = len(II)
    for isgm in range(nsgm):
      nLeg    = np.append(nLeg, ileg)
      if isgm==0:
        if ileg==0:
          LegNorm  = LegVnorm.copy()
      else:
        LegNorm = np.append(LegNorm, LegVnorm, axis=0)

# Compute half 1 of the segment:
      if isgm == 0:
        ii1 = II[isgm]
        ii2 = II[isgm]
        jj1 = JJ[isgm]
        jj2 = JJ[isgm]     
        L1  = 0.

      else:
        ii1 = II[isgm-1]
        ii2 = II[isgm]
        jj1 = JJ[isgm-1]
        jj2 = JJ[isgm]     
        if ii1 == ii2:
          sgm_north = True
        else:
          sgm_north = False

        if sgm_north:
          L1 = 0.5*DY[jj1,ii1]
        else:
          L1 = 0.5*DX[jj1,ii1]

# Find angle between the half-segment and the main section:
      sgmV     = np.array([ii2-ii1,jj2-jj1,0])
      sgmV_mgn = np.sqrt(np.dot(sgmV, sgmV))
      if sgmV_mgn < 1.e-30:
        Tht1 = 0.
      else:
        cosTht   = np.dot(LegV,sgmV)/(LegV_mgn * sgmV_mgn)
  #      print(f"isgm={isgm} sgmV_mgn={sgmV_mgn} cos={cosTht}")
        Tht1     = np.arccos(cosTht)
        RotSgm1  = np.sign(np.cross(sgmV, LegV)[2]) # direction of rotation: 
                                                    # segment to main section 
        Tht1     = RotSgm1*abs(Tht1)
     
# Find norm for half-segment 1:
      nrm_v1 = find_norm([ii1,jj1],[ii2,jj2], curve_ornt)
      nrm_v1 = np.expand_dims(nrm_v1, axis=0)

# Compute half 2:
      if isgm == nsgm-1:
        ii1 = II[isgm]
        ii2 = II[isgm]
        jj1 = JJ[isgm]
        jj2 = JJ[isgm]
        L2  = 0.
      else:
        ii1 = II[isgm]
        ii2 = II[isgm+1]
        jj1 = JJ[isgm]
        jj2 = JJ[isgm+1]
        if ii1 == ii2:
          sgm_north = True
        else:
          sgm_north = False

        if sgm_north:
          L2 = 0.5*DY[jj2,ii2]
        else:
          L2 = 0.5*DX[jj2,ii2]
        
# Find angle between the half-segment and the main section:
      sgmV     = np.array([ii2-ii1,jj2-jj1,0])
      sgmV_mgn = np.sqrt(np.dot(sgmV, sgmV))
      if sgmV_mgn < 1.e-30:
        Tht2 = 0.
      else:
        cosTht   = np.dot(LegV,sgmV)/(LegV_mgn * sgmV_mgn)
        Tht2     = np.arccos(cosTht)
        RotSgm2  = np.sign(np.cross(sgmV, LegV)[2]) # direction of rotation: 
                                                    # segment to main section 
        Tht2     = RotSgm2*abs(Tht2)

# Find norm for half-segment 2:
      nrm_v2 = find_norm([ii1,jj1],[ii2,jj2], curve_ornt)
      nrm_v2 = np.expand_dims(nrm_v2, axis=0)

      if len(vnrm1) == 0:
        vnrm1 = nrm_v1.copy()
        vnrm2 = nrm_v2.copy()
      else: 
        vnrm1  = np.append(vnrm1, nrm_v1, axis=0)
        vnrm2  = np.append(vnrm2, nrm_v2, axis=0)

      Lsgm1    = np.append(Lsgm1, L1)
      Lsgm2    = np.append(Lsgm2, L2)
#      Tht2Sct1 = np.append(Tht2Sct1, Tht1)
#      Tht2Sct2 = np.append(Tht2Sct2, Tht2)


# Combine legs at the pivoting points
#
#  At the corner (pivoting) points I[izr] = I[izr+1]
#      and J[izr]=J[izr+1]
#
#                  : 
#                  : Lsgm2[izr] = 0
#                  :
#           .......*---- nLeg[izr+1]=Leg2
# Lsgm1[izr+1]=0   | 
#                  | 
#       nLeg[izr] = Leg1
#
  nsgm = len(Isgm)
  Dij  = np.sqrt((np.diff(Isgm))**2 + (np.diff(Jsgm))**2)
  Idbl = np.where(Dij==0)[0]
  IZR1 = np.array([])
  IZR2 = np.array([])
  ccr  = -1
  for nI in range(len(Idbl)):
    ccr += 1
    izr = Idbl[nI]-ccr  # adjust for deleted rows
# Check that this is a pivoting point or a duplicate:
#  remove duplicated point, although there should not be
# duplicates except for pivots connecting the legs 
    if nLeg[izr] == nLeg[izr+1]: 
      print(f"define_segm: segm {izr} duplicated point, not a pivoting point")
# TODO: if there are duplicates - remove everywhere
    else:
      Isgm = np.delete(Isgm, izr)
      Jsgm = np.delete(Jsgm, izr)
      if Lsgm1[izr] == 0:
        izr1 = izr
        izr2 = izr+1
      else:
        izr1 = izr+1
        izr2 = izr
# Check: both half Lsgm should be 0 for adjacent segments 
      if Lsgm1[izr1] < 1.e-30:
        Lsgm1 = np.delete(Lsgm1, izr1)
      else:
        raise Exception("Pivot point: Lsgm1({izr1}) not 0")

      if Lsgm2[izr2] < 1.e-30:
        Lsgm2 = np.delete(Lsgm2, izr2)
      else:
        raise Exception("Pivot point: Lsgm2({izr2}) not 0")
     
# Designate pivoting points (connectors) as half of
# the two legs Numbers
      nLeg = np.delete(nLeg, izr+1)
      nLeg[izr] = 0.5*(nLeg[izr]+nLeg[izr+1])

      nrm1  = vnrm1[izr1]
      nrm2  = vnrm2[izr2]
      lnrm1 = np.sqrt(np.dot(nrm1,nrm1.transpose()))
      lnrm2 = np.sqrt(np.dot(nrm2,nrm2.transpose()))
      if lnrm1 < 1.e-30:
        vnrm1 = np.delete(vnrm1, izr1, axis=0)
      else:
        raise Exception("vnrm1({izr1}) not 0") 

      if lnrm2 < 1.e-30:
        vnrm2 = np.delete(vnrm2, izr2, axis=0)
      else:
        raise Exception("vnrm2({izr2}) not 0") 

  Isgm = Isgm.astype(int) 
  Jsgm = Jsgm.astype(int) 
  SGMT = SEGMENTS(Isgm, Jsgm, Lsgm1, Lsgm2, vnrm1, vnrm2, nLeg, \
                  curve_ornt, LegNorm)
#  SGMT = np.dtype({'curve_ornt': curve_ornt,
#                   'I_indx': Isgm,
#                   'J_indx': Jsgm,
#                   'half_Lsgm1': Lsgm1,
#                   'half_Lsgm2': Lsgm2,
#                   'half_norm1': vnrm1,
#                   'half_norm2': vnrm2})

  return SGMT

def find_norm(IJ1, IJ2, curve_ornt):
  """
    With respect to IJ1, IJ2 direction of a segment
    Use cross-product Vz(0,0,1) x V(i,j,0) to get 
    positive norm
    curve_ornt = pos: norm is to the left
               = neg: norm is to the right
  """
  Vz=np.array([0,0,1])
  V=np.array([IJ2[0]-IJ1[0], IJ2[1]-IJ1[1], 0])
# Check for 0 vector:
  V_len = np.sqrt(V[0]**2 + V[1]**2 + V[2]**2)
  if V_len < 1.e-32:
#    raise Exception(f"0 vector segment ({IJ1[0]} {IJ1[1]}), ({IJ2[0]} {IJ2[1]})")
    vnrm=np.array([0,0])
    return vnrm

  Unrm = np.cross(Vz,V)
  Unrm_len = np.sqrt(Unrm[0]**2 + Unrm[1]**2 + Unrm[2]**2)
  Unrm = Unrm/Unrm_len   # should be 0 or 1
  if curve_ornt[0:3] == 'neg':
    Unrm = -Unrm

  vnrm = Unrm[0:2] 
  vnrm = np.where(abs(vnrm) < 1.e-30, 0, vnrm)

  return vnrm

def match_indices(IIr,JJr,LONr,LATr,LON,LAT):
  """
    Match indices between 2 grids
    Given input indices IIr, JJr and LONr/LATr grid
    Find corresponding (closest) on LON/LAT grid
    LON, LAT - 1D (for Mercator grid) or 2D arrays
  """
  if type(IIr) == list:
    IIr = np.array(IIr)
  if type(JJr) == list:
    JJr = np.array(JJr)

  if len(IIr.shape) > 1 or len(JJr.shape) > 1:
    raise Exception("Dim of Input indices IIr/JJr should be <= 1 ")

  if len(LONr.shape) == 1:
    LATr, LONr = np.meshgrid(LATr, LONr, indexing='ij')
  elif len(LONr.shape) > 2:
    raise Exception("Dim of LONR/LATR should be 1 or 2")

  if len(LON.shape) == 1:
    LAT, LON = np.meshgrid(LAT, LON, indexing='ij') 
  elif len(LON.shape) > 2:
    raise Exception("Dim of LON/LAT should be 1 or 2")

  npnts = len(IIr)
  II = (np.zeros((npnts))-999).astype(int)
  JJ = (np.zeros((npnts))-999).astype(int)
  for ii in range(npnts):
    i0 = IIr[ii]
    j0 = JJr[ii]
    x0 = LONr[j0,i0]
    y0 = LATr[j0,i0]

    DST  = dist_sphcrd(y0, x0, LAT, LON)    
    jmin, imin = np.where(DST == np.min(DST))
    jmin = jmin[0]
    imin = imin[0]

    II[ii] = int(imin)
    JJ[ii] = int(jmin)

  return II, JJ

def find_closest_indx(LONr, LATr, LON, LAT, f_progr=False):
  """
    Given lon/lat of points in LONr, LATr 1D arrays
    Find corresponding (closest) on LON/LAT grid
    LON, LAT - 1D (for Mercator grid) or 2D arrays
  """
  if type(LONr) == list:
    LONr = np.array(LONr)
  if type(LATr) == list:
    LATr = np.array(LATr)

#  L2d = True
  if len(LON.shape) == 1:
    LAT, LON = np.meshgrid(LAT, LON, indexing='ij') 
#    L2d = False
  elif len(LON.shape) > 2:
    raise Exception("Dim of LON/LAT should be 1 or 2")

  npnts = len(LONr)
  II = (np.zeros((npnts))-999).astype(int)
  JJ = (np.zeros((npnts))-999).astype(int)
  for ii in range(npnts):
    if f_progr:
      if ii%50 == 0:
        fproc = float(ii)/npnts*100.
        print(f' processed {fproc:.1f}% ...')

    x0 = LONr[ii]
    y0 = LATr[ii]

    DST  = dist_sphcrd(y0, x0, LAT, LON)    
    jmin, imin = np.where(DST == np.min(DST))
    jmin = jmin[0]
    imin = imin[0]

    II[ii] = int(imin)
    JJ[ii] = int(jmin)

# Eliminate repitions, identical indices:
  IIn = list([II[0]])
  JJn = list([JJ[0]])
  for ii in range(1,len(II)):
    dd = np.sqrt( (II[ii]-IIn[-1])**2 + (JJ[ii]-JJn[-1])**2 )
    if dd < 1.e-6:
      continue
    IIn.append(II[ii])
    JJn.append(JJ[ii])

  II = np.array(IIn)
  JJ = np.array(JJn)

  return II, JJ

def orientation(A,B,C):
  """
    Given 3 points: A(ax,ay), B(bx,by), C(cx,cy)
    A,B,C - lists
   Given vector ab by two points: A(ax,ay) - start point of vector
    and B(bx,by) - end point, such that positive vector in X means
   bx-ax > 0
   and a point C(cx,cy)
   Does C lie on, to the left, to the right of ab?
   To find this out, find determinant of 
   | ax  ay  1 |
   | bx  by  1 | = | ax-cx  ay-cy |
   | cx  cy  1 |   | bx-cx  by-cy |
  
   From website of J. Shewchuk 
  
   Output: determinate D
   If point ay < by, then D<0 means point C is to the right
                          D>0 - point C is to the left
                          D=0 - point C is on the line AB
  """
  ax, ay = A
  bx, by = B
  cx, cy = C
  
  D = (ax-cx)*(by-cy)-(ay-cy)*(bx-cx);

  return D

def vector_projection(A, B):
  """
    Vector projection of A on vector B
    A, B are 2-element arrays
  """
  ndim = len(A.shape)
  if ndim == 1: A = np.expand_dims(A, axis=1)
  ndim = len(B.shape)
  if ndim == 1: B = np.expand_dims(B, axis=1)

  lA = np.sqrt(np.dot(A.transpose(),A))[0][0]
  lB = np.sqrt(np.dot(B.transpose(),B))[0][0]

  cosT = np.dot(A.transpose(), B)[0][0]/(lA*lB)
  prA  = lA*cosT
  tht_deg = np.arccos(cosT)*180./np.pi

  return prA, cosT, tht_deg

def construct_vector(pnt1, pnt2, vect21=True):
  """
    Given 2 pnts, construct a vector from (x1,y1) --> (x2,y2)
    pnt1, pnt2 - lists or np arrays (2 elements)
    vect21 = return a 2,1 array for linear algebra operations
    otherwise (2,0) array
  """
  x1, y1 = pnt1
  x2, y2 = pnt2
  Av = np.array([x2-x1,y2-y1])
  if vect21:
    Av = np.expand_dims(Av, axis=1)

  return Av

def shuffle3D_lon180(A3d, lon, conv360=True):
  """
    Global 3D array on Mercator grid, lon = 1D array
    lon = [-180, 180]
    want to rearrange to lon = [0, 360] such that -180/180 is in the middle of the data set
  """  
  if np.min(lon>=0.): 
    lon = np.where(lon>180, lon-360., lon)

  lon1 = lon[0]
  lon2 = lon[-1]
  if lon1 > lon2:
    raise Exception(f"Input lon needs to be increasing but lon1={lon1} lon2={lon2}")

  ilon0 = np.where(lon>=0)[0][0] 
  lonR  = np.concatenate((lon[ilon0:],lon[:ilon0]))
  if conv360: 
    lonR = np.where(lonR<0, lonR+360., lonR)

  A3dR = np.concatenate((A3d[:,:,ilon0:], A3d[:,:,:ilon0]), axis=2)

  return A3dR, lonR

def shuffle1D_lon180_to0360(lon, conv360=True):
  """
    Mercator grid, lon = 1D array
    lon = [-180, 180] with 0 somewhere in the middle of the grid
    i.e. the grid starts in the Western h/sphere and goes
    over the 0 meridian to the East
    want to rearrange to lon = [0, 360] such that discontinuity 
    -180/180 is in the middle of the grid
  """  
  if np.min(lon>=0.): 
    lon = np.where(lon>180, lon-360., lon)

  lon1 = lon[0]
  lon2 = lon[-1]
  if lon1 > lon2:
    raise Exception(f"Input lon has to be increasing but lon1={lon1} lon2={lon2}")

  ilon0 = np.where(lon>=0)[0][0] 
  lonR  = np.concatenate((lon[ilon0:],lon[:ilon0]))
  if conv360: 
    lonR = np.where(lonR<0, lonR+360., lonR)

  return lonR

def convert_polarXY_lonlat(X, Y, RE=6378137.0, E=0.08181919, SLAT=70., North=True, \
                           LON0_dir = -45., units='m'):
  """
  Transform polar X, Y coordinates (m) to geogr. lon/lat
  in N. hemisphere  and S. hemisphere (check signs for longitudes in S. hemisphere)
  The code is based on fortran version :
  https://github.com/nsidc/polarstereo-latlon-convert-fortran/blob/main/locate/locate.for
  as well as matlab version:
  https://www.mathworks.com/matlabcentral/fileexchange/32907-polar-stereographic-coordinate-transformation-map-to-lat-lon

  With modifications to make it work for ice conc grid 
     RE      : Radius of ellipsoid, WGS84
     E       : Eccentricity, WGS84
     E2      : Eccentricity squared
     SLAT    : Standard latitude for the SSM/I grids is 70 degrees
               Latitude of 0 distortion, i.e. where a projection plane 
               tangent to the Earth surface
     LON0_dir : orientation of the 0 longitude wth to X axis on polar grid
                sign = math sense (+ c/clockwise)

  Test:
  X=-3850e3 Y=5850e3 =>  lon = 168.35, lat=30.98  - left upper corner of NSDIC grid
  X=3750e3  Y=5850e3 =>  lat=31.37 lon=102.34
  X=3750e3  Y=-5350e3 => lat=34.35 lon=350.03

  """

# If the standard latitude is in S. Hemisph - switch signs
  sgn = 1.
  if SLAT < 0:
    sgn  = -1.
    SLAT = -SLAT
    LON0_dir = -LON0_dir
    X = -X # 
    Y = -Y #

# Keep everything in meters, input X, Y - meters
  cff = 1. 
  if not units == 'm':
    cff = 1.e3

  X = X*cff
  Y = Y*cff

# See dimensions:
  if isinstance(X, np.ndarray):
    ndim = len(X.shape)
  else:
    ndim = 0

  deg2rad = np.pi/180.
  E2    = E**2
  phi0  = SLAT*deg2rad
  lmbd0 = LON0_dir*deg2rad

  rho  = np.sqrt(X**2 + Y**2)

  if ndim == 0:
    if rho < 0.1:
      alat = sgn*90.
      alon = 0.
      return alon, alat

  cm  = np.cos(phi0)/(np.sqrt(1.-E2*((np.sin(phi0))**2)))
  a0  = np.tan((np.pi/4.) - (phi0/2.))
  a1  = ((1.-E*np.sin(phi0))/(1.0+E*np.sin(phi0)))**(E/2.)
  t_c = a0/a1 

  if (abs(SLAT-90.) < 1.e-5): 
    TT = rho*np.sqrt((1.+E)**(1.+E)*(1.-E)**(1.-E))/2./RE
  else:
    TT = rho*t_c/(RE*cm)

# Find phi using a series:
  chi  = np.pi/2. - 2.*np.arctan(TT)
  b1   = (E2/2. + 5.*(E2**2)/24. + (E**6)/12 + 13.*(E**8)/360)*np.sin(2.*chi)
  b2   = (7*(E**4)/48. + 29.*(E**6)/240. + 811*(E**8)/11520.)*np.sin(4.*chi)
  b3   = (7.*(E**6)/120.+81.*(E**8)/1120.)*np.sin(6.*chi)
  b4   = (4279.*(E**8)/161280.)*np.sin(8.*chi)
  phi  = chi + b1 + b2 + b3 + b4

# Adjust sign and long rotation
  phi  = sgn*phi
  lmbd = np.arctan2(Y, X) - sgn*lmbd0

  lon = lmbd*180./np.pi
  lat = phi*180./np.pi

# Make long -180, 180
  lon = ((lon+180.)%360.) - 180.

  return lon, lat


def connect_segments(Isgm, Jsgm):
  """
    Connect N segments together
    Isgm = list of N segments i-coord
    Jsgm = list of N segments j-coord

    Segments orientation is checked to have continuous transition from segm 1 to segm 2
  """
  nsegm = len(Isgm)

  II = Isgm[0]
  JJ = Jsgm[0]

  if isinstance(II, list):
    II = np.array(II)
    JJ = np.array(JJ)

  for isgm in range(1,nsegm):
    I1 = Isgm[isgm]
    J1 = Jsgm[isgm]

    if isinstance(I1, list):
      I1 = np.array(I1)
      J1 = np.array(J1)

    d1 = np.sqrt((I1-II[-1])**2 + (J1-JJ[-1])**2)
    if d1[0] > d1[-1]:
      I1 = np.flip(I1)
      J1 = np.flip(J1)

    II = np.append(II, I1)
    JJ = np.append(JJ, J1)

  return II, JJ


def convert_nptime_to_datenum(Time):
  """
    Convert an array of numpy datetime64 to date numbers
    datetime data are typically derived from netcdf Time variables
  """
  import datetime
  from datetime import datetime
  import mod_time as mtime

  nrec = len(Time)
  TM = np.zeros((nrec))
  for irc in range(nrec):
    dmm = Time[irc].astype('datetime64[s]')
    dmm = dmm.astype(datetime)
    yr  = dmm.year
    month = dmm.month
    mday  = dmm.day
    hour  = dmm.hour
    mint  = dmm.minute
    TM[irc] = mtime.datenum([yr,month,mday,hour,mint])

  return TM









