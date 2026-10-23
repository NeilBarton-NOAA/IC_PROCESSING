"""
  Sea water state functions 
    based on UNESCO algorithms
    The code is based on  Phil Morgan 1993, CSIRO
    Matlab version


  2026: 
  added Ocean Freezing T similar to MOM6
  see: MOM6/src/equation_of_state/MOM_TFreeze.F90

  Dmitry Dukhovskoy, NOAA NWS EMC
"""
import numpy as np
import sys
#import pdb

def sw_press(zz,LAT):
  """
  Compute pressure (dBar, Pa) for given depths in m
  LAT - local latitudes (degrees)

  REFERENCE: Saunders, P.M., 1981
  "Practical conversion of Pressure to Depth"
  JPO, 11, 573-574

  The code is similar to the sw_pres.m developed by
  Phil Morgan 1993, CSIRO

  zz - scalar, or arrays
  """

#  if np.isscalar(zz):
#    lat0 = LAT
#  elif isinstance(zz,np.ndarray):
#    if np.dim(zz)==1:
#      lat0=LAT
#    else:
#      nrws = zz.shape[0]
#      nclm = zz.shape[1]
#      
#    for ii in range(nclm):
#      lat0 = LAT[:,ii]
  zz=np.abs(zz)  
##      breakpoint()
# Check shapes 
  zdim = 0
  if isinstance(zz,np.ndarray):
    if not isinstance(LAT,np.ndarray):
      print('ERROR: Depth Z is array and LAT is not array')
      raise Exception('Inputs must be both arrays or scalars')
    else:
      zdim = zz.ndim
      ldim = LAT.ndim
      if zdim != ldim:
        print('ERROR: Z dim={0}, Lat dim={1}, should be equal'.format(zdim,ldim))
        raise Exception('dimensions disagree Z and Lat')

  deg2rad = np.pi/180.
  X    = np.sin(np.abs(LAT)*deg2rad)
  C1   = 5.92e-3 + (X**2)*5.25e-3
  pres_db = ((1.-C1)-np.sqrt(((1.-C1)**2)-(zz*8.84e-6)))/(4.42e-6) # dbar
 
  # Make sure pressure = 0 at the surface
  if zdim > 0:
    pres_db = np.where(zz < 1.e-5, 0., pres_db) 

  pres_pa = pres_db*1.e4  # Pa

  return pres_db, pres_pa

def adiab_Tgrad(S,T,P):
  """
    Calculates adiabatic T gradient as per UNESCO 1983
    all input fields must have same dim

    Possible types: scalar (float) or 2D numpy arrays
    
    S  = salinity, psu PSS-78
    T  = temperature, C dgr ITS-90
    P  = pressure, db

    Based on P. Morgan code
  """
  
# Check dim:
  if isinstance(S, np.ndarray):
#    print('type S={0}'.format(type(S))
    ndim = S.ndim
    jS  =  S.shape[0]
    jT = T.shape[0]
    jP = P.shape[0]

    iS = 0
    iT = 0
    iP = 0
    if ndim==2:
      iS = S.shape[1]
      iT = T.shape[1]
      iP = P.shape[1]

    if jS != jT or jS != jP or iS != iP or iS != iT:
      print('adiab_Tgrad: ERROR dimensions of T,S,P do not match')
      print('{0}x{1} {2}x{3} {4}x{5}'.format(jS,iS,jT,iT,jP,iP))
      raise Exception ('error dim ')

  T68 = 1.00024*T

  a0 =  3.5803e-5
  a1 =  8.5258e-6
  a2 = -6.836e-8
  a3 =  6.6228e-10

  b0 =  1.8932e-6
  b1 = -4.2393e-8
  
  c0 =  1.8741e-8
  c1 = -6.7795e-10
  c2 =  8.733e-12
  c3 = -5.4481e-14

  d0 = -1.1351e-10
  d1 =  2.7759e-12

  e0 = -4.6206e-13
  e1 =  1.8676e-14
  e2 = -2.1687e-16

  adtg =   a0 + (a1 + (a2 + a3*T68)*T68)*T68 \
         + (b0 + b1*T68)*(S-35.) \
         + ((c0 + (c1 + (c2 + c3*T68)*T68)*T68) \
         + (d0 + d1*T68)*(S-35.))*P \
         + (e0 + (e1 + e2*T68)*T68)*P*P

  return adtg


#
def sw_ptmp(S,T,P1,Pref):
  """
    Potential T at reference pressure Pref adiabatically brought
    from Pressure=P1
    Based on UNESCO 1983 report - Phil Morgan Matlab code
    all variables have same dimensions (e.g. 2D arrays)
    S = salinity psu PSS-78
    T = temperature C deg ITS-90
    P1 = pressure dbar
    Pref = pressure reference dbar
  """

  if np.min(Pref)<0. or np.min(P1)<0:
    print('ERROR: Pressure must be >0 min Pref={0} P1={1}'.\
          format(np.min(Pref),np.min(P1)))
    raise Exception(' Pressure error') 

# theta 1
  del_P   = Pref - P1
  adtg    = adiab_Tgrad(S,T,P1)
  del_th  = del_P*adtg
  th      = T*1.00024 + 0.5*del_th
  q       = del_th

# theta 2
  adtg = adiab_Tgrad(S,th/1.00024,P1+0.5*del_P)
  del_th = del_P*adtg
  th     = th + (1.-1./np.sqrt(2.))*(del_th-q)
  q      = (2. - np.sqrt(2.))*del_th + (-2.+3./np.sqrt(2.))*q

# theta 3
  adtg   = adiab_Tgrad(S,th/1.00024,P1+0.5*del_P)
  del_th = del_P*adtg
  th     = th + (1.+1./np.sqrt(2.))*(del_th-q)
  q      = (2.+np.sqrt(2.))*del_th + (-2.-3./np.sqrt(2.))*q

# theta 4
  adtg   = adiab_Tgrad(S,th/1.00024,P1+del_P)
  del_th = del_P*adtg
  ptmp   = (th + (del_th - 2.*q)/6.)/1.00024

  return ptmp

#
#
def sw_smow(T):
  """
    Density of standard mean ocean water (pure water)
    UNESCO EOS 1980    

    Input T C deg ITS-90

    Output rho_S0 - kg/m3
  """
  a0 =  999.842594
  a1 =  6.793952e-2
  a2 = -9.095290e-3
  a3 =  1.001685e-4
  a4 = -1.120083e-6
  a5 =  6.536332e-9

  T68 = T*1.00024
  rho_S0 = a0 + (a1 + (a2 + (a3 + (a4 + a5*T68)*T68)*T68)*T68)*T68

  return rho_S0



def sw_dens0(S,T):
  """
    Sea water density at pressure=0, i.e. atmospheric pressure
    UNESCO 1983 (EOS 1980) polynomial

    S = salinity psu PSS-78
    T = temp. C deg ITS-90

    output: 
      rhoP0 kg/m3
    code is based on Phil Morgan's Matlab code

  """
  T68 = T*1.00024

# UNESCO 1983 eqn(13) p17
  b0 =  8.24493e-1
  b1 = -4.0899e-3
  b2 =  7.6438e-5
  b3 = -8.2467e-7
  b4 =  5.3875e-9

  c0 = -5.72466e-3
  c1 =  1.0227e-4
  c2 = -1.6546e-6

  d0 =  4.8314e-4

  rho_P0 = sw_smow(T) + (b0 + (b1 + (b2 + (b3 + b4*T68)*T68)*T68)*T68)*S \
                    + (c0 + (c1 + c2*T68)*T68)*S*np.sqrt(S) + d0*S**2

  return rho_P0


def sw_seck(S,T,P):
  """
    Secant bulk modulus (K) of sea water using equation of state 1980
    UNESCO polynomial

    S = salinity PSS-78
    T = temp  ITS-90
    P = pressure, db

    output
    K = Secant bulk modulus, bars
  """

#--------------------------------------------------------------------
#
# COMPUTE COMPRESSION TERMS
#
#--------------------------------------------------------------------

  P = P/10.  #convert from db to atmospheric pressure units
  T68 = T * 1.00024

# Pure water terms of the secant bulk modulus at atmos pressure.
# UNESCO eqn 19 p 18

  h3 = -5.77905E-7
  h2 = +1.16092E-4
  h1 = +1.43713E-3
  h0 = +3.239908   #[-0.1194975];

  AW  = h0 + (h1 + (h2 + h3*T68)*T68)*T68

  k2 =  5.2787E-8
  k1 = -6.12293E-6
  k0 =  8.50935E-5   #[+3.47718E-5];

  BW  = k0 + (k1 + k2*T68)*T68

  e4 = -5.155288E-5
  e3 =  1.360477E-2
  e2 = -2.327105
  e1 =  148.4206
  e0 =  19652.21    #[-1930.06];

  KW  = e0 + (e1 + (e2 + (e3 + e4*T68)*T68)*T68)*T68   #eqn 19

#--------------------------------------------------------------------
#
# SEA WATER TERMS OF SECANT BULK MODULUS AT ATMOS PRESSURE.
#
#--------------------------------------------------------------------

  j0 = 1.91075E-4 

  i2 = -1.6078E-6
  i1 = -1.0981E-5
  i0 =  2.2838E-3

  SR = np.sqrt(S)

  A  = AW + (i0 + (i1 + i2*T68)*T68 + j0*SR)*S

  m2 =  9.1697E-10
  m1 = +2.0816E-8
  m0 = -9.9348E-7

  B = BW + (m0 + (m1 + m2*T68)*T68)*S   # eqn 18

  f3 =  -6.1670E-5
  f2 =  +1.09987E-2
  f1 =  -0.603459
  f0 = +54.6746

  g2 = -5.3009E-4
  g1 = +1.6483E-2
  g0 = +7.944E-2

  K0 = KW + (f0 + (f1 + (f2 + f3*T68)*T68)*T68 \
        +   (g0 + (g1 + g2*T68)*T68)*SR)*S      # eqn 16
  K = K0 + (A + B*P)*P  # eqn 15

  return K


# Calculate sea water density at reference pressure 
def sw_dens(S,T,P):
  """
    Computes seawater density using UNESCO 1983 (EOS 80) polynomial

    All input fields must have same dimensions (or scalar)
    S = salinity PSS-78
    T = temp C deg, ITS-90
    P = pressure (>0 !) db

    output density rho kg/m3
  """

  rhoP0 = sw_dens0(S,T)
  K     = sw_seck(S,T,P)
  P     = P/10.   # db --> bar
  rho   = rhoP0/(1.-P/K)

  return rho

# Ocean freezing T:
def freezing_temp_linear(Sal, pres):
  """
  Compute freezing point potential temperature of seawater.

  Parameters:
      S (float): Salinity [PSU]
      pres (float): Pressure [Pa]

  Returns:
      float: Freezing temperature [degC]

  From MOM6:
!> This subroutine computes the freezing point potential temperature
!! [degC] from salinity [ppt], and pressure [Pa] using the expression
!! from Millero (1978) (and in appendix A of Gill 1982), but with the of the
!! pressure dependence changed from 7.53e-8 to 7.75e-8 to make this an
!! expression for potential temperature (not in situ temperature), using a
!! value that is correct at the freezing point at 35 PSU and 5e6 Pa (500 dbar).

  """
  # Constants
  cS1     = -0.0575        # A term in the freezing point fit [degC PSU-1]
  cS3_2   = 1.710523e-3    # A term in the freezing point fit [degC PSU-3/2]
  cS2     = -2.154996e-4   # A term in the freezing point fit [degC PSU-2]
  dTFr_dp = -7.75e-8       # Derivative of freezing point with pressure [degC Pa-1]

  # Convert to numpy arrays (handles scalar or array input)
  Sal = np.asarray(Sal)
  pres = np.asarray(pres)

  # Ensure non-negative salinity inside sqrt
  S_sqrt = np.sqrt(np.maximum(Sal, 0.0))

  # Compute freezing temperature
  Tfrz = Sal * (cS1 + (cS3_2 * S_sqrt + cS2 * Sal)) + dTFr_dp * pres

  return Tfrz


