"""
  Functions/subroutines for bilinear interpolation
  Dmitry Dukhovskoy NOAA NWS OMD
"""
import numpy as np

def basisFn_RectRef():
  """
    Define basis functions on a reference rectangle (square)
    (-1,-1), (1,-1) (1,1) (-1,1)

    For each basis fn phi1, ..., phi4 solve a system of linear eqns
    such that phi1 = 1 at (x1,y1) and 0 at the other vertices, etc.
    a + b*x1 + c*y1 + d*x1*y1 = 1
    a + b*x2 + c*y2 + d*x2*y2 = 0
    etc

    M*A = X, M is [1 -1 -1 1; 
                   1 1 -1 -1; etc, 
    A is (a; b; c; d), X=[1, 0, 0, 0] for phi1 
    See p. 83-90, Gockenbach, "Understanding Finite-Element Methods" textbook
    inverse of M is known
  """
  Minv = 1./4.*np.array([[1, 1, 1, 1],[-1,1,1,-1],[-1,-1,1,1],[1,-1,1,-1]])
  X1 = np.transpose(np.array([[1,0,0,0]]))
  X2 = np.transpose(np.array([[0,1,0,0]]))
  X3 = np.transpose(np.array([[0,0,1,0]]))
  X4 = np.transpose(np.array([[0,0,0,1]]))

  phi1 = np.dot(Minv,X1).squeeze()
  phi2 = np.dot(Minv,X2).squeeze()
  phi3 = np.dot(Minv,X3).squeeze()
  phi4 = np.dot(Minv,X4).squeeze()

  return phi1, phi2, phi3, phi4

def sort_gridcell_indx(II, JJ, fsens='positive'):
  """
    Alternative: mmisc.reorder_polygon(xv, yv, indx=False)

    In case when rectangular indices are in random order
    need to put them in c/clckwize (positive) or clockwize (negative)
    order

    Note: works only for simple cases for quadrilateral recatngular box !!!
    i.e. there are only 2 different pairs of values for I and J
    this does not work for "trapezoid" shaped box and similar

    assumed: indices represent vertices of a grid box:
    i = 0,1,2,3

     II[i]JJ[i]     II[..]JJ[..]
     *----------------*
     |                |
     |                |
     |                |
     *----------------*
     II[..]JJ[..]     II[..]JJ[..]

   e.g.: II = [10,11,10,11] J=[23, 22, 22, 23]
     want:  II=[10,11,11,10] J=[22,22,23,23]
  """
  if isinstance(II, list):
    II = np.array(II)
  if isinstance(JJ, list):
    JJ = np.array(JJ)
    
# Start from smallest i,j:
  D    = np.sqrt(II**2 + JJ**2)
  isrt = np.argsort(D)
  I0   = II.copy()
  J0   = JJ.copy()  

  ix1         = I0[isrt[0]]
  jx1         = J0[isrt[0]]
  I0[isrt[0]] = -999
  J0[isrt[0]] = -999

  if fsens == 'positive':
# c/clockwise sense:
    imm     = np.argwhere(J0 == jx1)[0][0]
    isrt[1] = imm
    ix2     = I0[imm]
    jx2     = J0[imm]
    I0[imm] = -999
    J0[imm] = -999

    imm     = np.argwhere(I0 == ix2)[0][0]
    isrt[2] = imm
    ix3     = I0[imm]
    jx3     = J0[imm]
    I0[imm] = -999
    J0[imm] = -999
 
    imm     = np.argwhere(J0 == jx3)[0][0]
    isrt[3] = imm

  else:
# clockwise sense:
    imm     = np.argwhere(I0 == ix1)[0][0]
    isrt[1] = imm
    ix2     = I0[imm]
    jx2     = J0[imm]
    I0[imm] = -999
    J0[imm] = -999

    imm     = np.argwhere(J0 == jx2)[0][0]
    isrt[2] = imm
    ix3     = I0[imm]
    jx3     = J0[imm]
    I0[imm] = -999
    J0[imm] = -999

    imm     = np.argwhere(I0 == ix3)[0][0]
    isrt[3] = imm

  IIsrt = II[isrt]
  JJsrt = JJ[isrt]

  return IIsrt, JJsrt

def phi_x0y0(phi,x1,y1):
  """
    Calculate basis fn phi at pnt(x1,y1)
  """
  bb1 = phi[0] + phi[1]*x1 + phi[2]*y1 + phi[3]*x1*y1

  return bb1


def map_x2xhat(XV,YV,x0,y0):
  """
    Map a point from a quadrilateral (x1,y1), ..., (x4,y4) 
    to a reference square
    (-1,-1), (1,-1), (1,1), (-1,1)
    The mapping is of the form:
    xhat = a1+a2*x+a3*y+a4*x*y = Alf*XY0
    yhat = b1+b2*x+b3*y+b4*x*y = Bet*XY0
    Solve 4 equations to find coefficients:
    a1+a2*x1+a3*y1+a4*x1*y1 = -1 
    a1+a2*x2+a3*y2+a4*x2*y2 = 1 
    ...

    AA*Alf = Cx

   same for b coefficients

   Input: 
     XV, YV - X, Y coordinates of 4 nodes, as 1D vectors of a grid cell
              encompassing the interpolation point
     x0,y0 - coordinates of the interpolation point (same unites as X,Y)
  """
  AA=np.array([[1, XV[0], YV[0], XV[0]*YV[0]],\
               [1, XV[1], YV[1], XV[1]*YV[1]],\
               [1, XV[2], YV[2], XV[2]*YV[2]],\
               [1, XV[3], YV[3], XV[3]*YV[3]]])

#  AA = np.float64(AA)
  AAinv = np.linalg.inv(AA)
  Cx = np.transpose(np.array([[-1.,1.,1.,-1.]]))
  Alf = np.dot(AAinv,Cx)

  Cy = np.transpose(np.array([[-1.,-1.,1.,1.,]]))
  Bet = np.dot(AAinv,Cy)


  XY0 = np.transpose(np.array([[1,x0,y0,x0*y0]]))
  xhat = np.dot(Alf.transpose(),XY0)[0][0]
  yhat = np.dot(np.transpose(Bet),XY0)[0][0]

  return xhat, yhat

def map_xhat2x(XV, YV, xht, yht):
  """
    Inversing mapping from a reference quadrilateral ---> original
    quadrilateral with vertices XV, YV
    The mapping is:
    a1+a2*xhat+a3*yhat+a4*xhat*yhat = x
    b1+b2*xhat+b3*yhat+b4*xhat*yhat = y

    Solve a system of lin eq:
    a1 + a2*s1 + a3*t1 + a4*s1*t1 = x1
    a1 + a2*s2 + a3*t2 + a4*s2*t2 = x2
    ...

    b1 + b2*s1 + b3*t1 + b4*s1*t1 = y1
    b1 + b2*s2 + b3*t2 + b4*s2*t2 = y2
    ...

  """
  XV    = np.expand_dims(XV, axis=0)
  YV    = np.expand_dims(YV, axis=0)

  AA = np.array([[1, -1, -1,  1],\
                 [1,  1, -1, -1],\
                 [1,  1,  1,  1],\
                 [1, -1,  1, -1]])

  AAinv = np.linalg.inv(AA)
  Alf   = np.matmul(AAinv,XV.transpose())
  Bet   = np.matmul(AAinv,YV.transpose())
  ST0   = np.transpose(np.array([[1, xht, yht, xht*yht]]))
  xx0   = np.dot(Alf.transpose(), ST0)[0][0]
  yy0   = np.dot(Bet.transpose(), ST0)[0][0]

  return xx0, yy0

def bilin_interp(phi1,phi2,phi3,phi4,xht,yht,HT):
  """
  Bilinear interpolation H(x,y) = sum(H(x1,y1)*psi_i(x,y)), i=1,...,4
  Input: 
    phi1,phi2,phi3,phi4 - linear basis functions
    xht, yht - interpolated point mapped onto a reference square
               use map_x2xhat   
  """
  p1x = phi_x0y0(phi1,xht,yht) # Phi1 at pnt xht,yht
  p2x = phi_x0y0(phi2,xht,yht)
  p3x = phi_x0y0(phi3,xht,yht)
  p4x = phi_x0y0(phi4,xht,yht)

  Hi = HT[0]*p1x + HT[1]*p2x + HT[2]*p3x + HT[3]*p4x

  return Hi

def bilin_interp1D(phi1,phi2,phi3,phi4,xht,yht,HT):
  """
  Bilinear interpolation H(x,y) = sum(H(x1,y1)*psi_i(x,y)), i=1,...,4
  For 1 D inpyt arrays combined in HT = [[1:n],[1:n],[1:n],[1:n]]
  where n - number of observations at each box vertex, e.g.
  n = number of vertical layers in a model   

  Input: 
    phi1,phi2,phi3,phi4 - linear basis functions
    xht, yht - interpolated point mapped onto a reference square
               use map_x2xhat   
  """
  p1x = phi_x0y0(phi1,xht,yht) # Phi1 at pnt xht,yht
  p2x = phi_x0y0(phi2,xht,yht)
  p3x = phi_x0y0(phi3,xht,yht)
  p4x = phi_x0y0(phi4,xht,yht)

  if not HT.ndim == 2:
    raise Exception('HT should be N x 4 array')

  npnt = HT.shape[1]
  nlv  = HT.shape[0]
  PX  = np.array([[p1x, p2x, p3x, p4x]]).transpose()
  Hi = np.matmul(HT,PX)
#  Hi = (HT[:,0]*p1x + HT[1]*p2x + HT[2]*p3x + HT[3]*p4x)[0]
  Hi = np.squeeze(Hi)

  return Hi

def blinrint(xht,yht,HT):
  """
  Bilinear interpolation H(x,y) = sum(H(x1,y1)*psi_i(x,y)), i=1,...,4
  Input: 
    xht, yht - interpolated point mapped onto a reference square
               use map_x2xhat   
    HT       - 1D array of node values
  """
# Find basis functions for a reference square:
  phi1,phi2,phi3,phi4 = basisFn_RectRef()

  p1x = phi_x0y0(phi1,xht,yht) # Phi1 at pnt xht,yht
  p2x = phi_x0y0(phi2,xht,yht)
  p3x = phi_x0y0(phi3,xht,yht)
  p4x = phi_x0y0(phi4,xht,yht)

  Hi = (HT[0]*p1x + HT[1]*p2x + HT[2]*p3x + HT[3]*p4x)

  return Hi

def rotate_box(XV, YV, x0c, y0c):
  """
    For tilted, rotated quadrilaterals mapping does not work well
    Try to rotate the box to align it with X/Y axes
    XV, YV coordinates for 4 vertices
    x0c, y0c - coordinates of the grid point inside the box where values are interpolated
               if this is not needed then can provide any x0c, y0c inside the box
    rotate around the 1st vertex
  """
  XR = XV.copy()
  YR = YV.copy()
#  x00c = x0c
#  y00c = y0c
  x0c = x0c-XV[0]
  y0c = y0c-YV[0]
  XV = XV-XV[0]
  YV = YV-YV[0]

  x1, x2, x3, x4 = XV
  y1, y2, y3, y4 = YV
  Av  = mmisc1.construct_vector([x1,y1],[x2, y2])
  Bv  = mmisc1.construct_vector([x1,y1],[x1+1.,y1]) # horizontal axis
  prA, cosT, tht = mmisc1.vector_projection(Av, Bv)
  if Av[1,0] < 0.: tht = -tht

  tht_rad = tht*np.pi/180.
  RR = np.array([[np.cos(tht_rad), np.sin(tht_rad)],[np.sin(-tht_rad), np.cos(tht_rad)]])

  for ii in range(4):
    U  = np.array([[XV[ii]],[YV[ii]]])
    Ur = np.matmul(RR,U).squeeze()
    XR[ii] = Ur[0]
    YR[ii] = Ur[1]

  U = np.array([[x0c],[y0c]])
  Ur = np.matmul(RR,U).squeeze()
  x0r = Ur[0]
  y0r = Ur[1]

  return XR, YR, x0r, y0r

def shift_longitudes(lons, ref_lon=360):
  """
  Normalize longitudes (degrees East) to a consistent range near a reference.
  If ref_lon is given, wrap longitudes so they are within +/-180dgr of ref_lon.
  The function shifts all longitudes so that they are around and close to the reference 
  longitue (ref_lon = x0).
  Default: normalize to [0, 360)
  For interpolation to avoid discontinuities near 0/360 or -180/180
  e.g. ref_lon = 185 and lons are [-180,180] --> the range
       does not include the ref_lon, after shifting the range will be [0,360] 
       to include the ref_lon
  """
  lons = np.array(lons, dtype=float)
  if ref_lon == 360:
    lons = np.mod(lons, 360.0)
  else:
    # Shift so that longitudes are within +/-180° of ref_lon
    lons = (lons - ref_lon + 180.0) % 360.0 - 180.0 + ref_lon
  return lons

def check_convex(X, Y):
  """
    Check if a quadrilateral is convex:
    no "batterfly" (self-intersecting)
    no concave quad
    Important for bilenear interpolation
  """
  cross = []
  for ii in range(len(X)):
    dx1 = X[(ii+1)%4] - X[ii]
    dy1 = Y[(ii+1)%4] - Y[ii]
    dx2 = X[(ii+2)%4] - X[(ii+1)%4]
    dy2 = Y[(ii+2)%4] - Y[(ii+1)%4]
    cross.append(dx1*dy2 - dy1*dx2)

  convex =  np.all(np.array(cross) > 0) or np.all(np.array(cross) < 0)

  return convex

def strictly_convex_quad(X, Y, eps_tol=1e-4):
  """
    More stric check of convexivity for a quadrilateral
    guarantees that none of the vertices
    are on the same line
    Scale the corss product to get the nearly-collinear cases 
    when the cross ~0 relative to large edge lengths
    For scaled cross for coliniear lines is ~1e-6 - 1e-7
  """
  signs = []
  for i in range(4):
    dx1 = X[(i+1)%4] - X[i]
    dy1 = Y[(i+1)%4] - Y[i]
    dx2 = X[(i+2)%4] - X[(i+1)%4]
    dy2 = Y[(i+2)%4] - Y[(i+1)%4]

    cross = dx1*dy2 - dy1*dx2

    # scale factor ~ product of edge lengths
    scale = (dx1*dx1 + dy1*dy1)**0.5 * (dx2*dx2 + dy2*dy2)**0.5

    if scale == 0:
      return False

    cross_scaled = cross / scale
    #print(f"cross={cross}  scaled={cross_scaled}")

    if abs(cross_scaled) < eps_tol:
      #print(f"Convex test: Rejected < {eps_tol}")
      return False  # catches if --*---*  

    signs.append(np.sign(cross_scaled))

  convex = all(s > 0 for s in signs) or all(s < 0 for s in signs)
 
  return convex 

def lonlat2xy_wrtX0(XX, YY, x0, y0):
  """
    Convert lon/lat coordinates of 4 quadrilateral vertices XX, YY
    into cartisian coordinates (meters) wrt to a ref point x0, y0
    so that x0,y0 --> 0,0

    Note: not recommended to use x0 y0 if it is outside the quadrilateral
    as it warps the quadrilateral 
    better to use centroid as a reference
    or some other reference pnt 
    and XX, YY ---> +/- distance from 0,0
  """
  import mod_misc1 as mmisc1

  # The method does not work near the poles!!!
  # Use lonlat2xy_enu instead
  assert(abs(y0) < 89.95), f"This method does not work near the poles"
  assert(np.max(abs(YY)) < 89.95), f"This method is not robust near the poles"

  # Normalize reference longitude
  xrf = ((x0 + 180) % 360) - 180
  yrf = y0

  # Define positive east and north vectors
  # originating at the reference point
  # Construct east/north direction vectors
  Ev = mmisc1.construct_vector([x0, y0], [x0+10., y0])
  Nv = mmisc1.construct_vector([x0, y0], [x0, y0+10.])
  XV = np.zeros(4)
  YV = np.zeros(4)

  for ipp in range(4):
    xx = ((XX[ipp] + 180) % 360) - 180  # normalize longitude
    yy = YY[ipp]

    # Vector from reference to point
    Av = mmisc1.construct_vector([x0, y0], [xx, yy])

    prAE, *_ = mmisc1.vector_projection(Av, Ev)
    prAN, *_ = mmisc1.vector_projection(Av, Nv)

    dlx = mmisc1.dist_sphcrd(yrf, xrf, yrf, xx)
    dly = mmisc1.dist_sphcrd(yrf, xrf, yy, xrf)

    XV[ipp] = np.sign(prAE) * dlx
    YV[ipp] = np.sign(prAN) * dly

  return XV, YV

from pyproj import Transformer

def make_lonlat2xy_transformer(x0, y0):
  """
    Local azimuthal equidistant projection centered at reference
    To convert lon/lat near the poles (or anywhere) into Cartesian distance 
    use a local tangent-plane coordinate system, also known as:
    East-North-Up (ENU)

    Create AEQD transformer centered at reference 
    To use for multiple transformations to spped up the process
  """
  from pyproj import Proj, Transformer
  return Transformer.from_crs(
    "EPSG:4326",
    f"+proj=aeqd +lat_0={y0} +lon_0={x0} +datum=WGS84",
    always_xy=True
  )


def lonlat2xy_enu(XX, YY, x0, y0):
  """
    Local azimuthal equidistant projection centered at reference
    To convert lon/lat near the poles (or anywhere) into Cartesian distance 
    use a local tangent-plane coordinate system, also known as:
    East-North-Up (ENU)

    Note: slow, do not use it inside for loops for many points
    instead: use   make_lonlat2xy_transformer to build the transformer outside the loop
             then use without rebuilding: xv, yv = transformer.transform(XX, YY)
  """
  from pyproj import Proj, Transformer

  proj_local = Proj(proj='aeqd', lat_0=y0, lon_0=x0, datum='WGS84')

  transformer = Transformer.from_proj("epsg:4326", proj_local, always_xy=True)

  x, y = transformer.transform(XX, YY)
  return x, y

def lonlat2xy_fast(xx, yy, lon0, lat0):
  """
    Approximate local ENU (meters) for small boxes (~50 km)
    Very fast alternative to pyproj Transformer for small regions.
  """
  R_earth = 6371000.0  # meters
  # Compute difference in longitude, wrapped to [-180, 180]
  dlon = (xx - lon0 + 180) % 360 - 180
  dlat = yy - lat0

  # Convert to local ENU distances
  x = dlon * np.cos(np.deg2rad(lat0)) * np.pi / 180 * R_earth
  y = dlat * np.pi / 180 * R_earth

  return x, y

def lonlat2xy_pnt(xx, yy, x0, y0):
  """
    Convert lon/lat coordinates of 1 point
    into cartisian coordinates (meters) wrt to a reference point (x0, y0)

    so that x0,y0 --> 0,0
    and XX, YY ---> +/- distance from 0,0

    Not robust near poles and lonh discontinuities
    use lolat2xy_enu
  """
  import mod_misc1 as mmisc1

  xrf = x0
  if xrf < -180.:
    xrf = xrf+360.
  yrf = y0

# Define positive east and north vectors
# originating at the reference point
  Ev  = mmisc1.construct_vector([x0, y0], [x0+10., y0])
  Nv  = mmisc1.construct_vector([x0, y0], [x0, y0+10.])
  Av  = mmisc1.construct_vector([x0, y0],[xx, yy])
  prAE, cosT, tht = mmisc1.vector_projection(Av, Ev)
  prAN, cosT, tht = mmisc1.vector_projection(Av, Nv) 
  dlx = mmisc1.dist_sphcrd(yrf, xrf, yrf, xx)
  dly = mmisc1.dist_sphcrd(yrf, xrf, yy, xrf)
  XV  = np.sign(prAE)*dlx
  YV  = np.sign(prAN)*dly

  return XV, YV


  
