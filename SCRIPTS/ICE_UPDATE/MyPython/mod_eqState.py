"""
  Eq. of state used in HYCOM
  check blkdat.input to see what eq. of state is used
  sigver=6  !17-term sigma-2
  sigver=5  !17-term sigma-0

  in hycom source, see header file
  stmt_fns.h
  
"""
import numpy as np

def sig2_17t_db(t,s,sigR=2000.):
  """
  17-term equation of state used in RTOFS, GOFS3.2, 3.5
  sigR - reference pressure / depth (db)
  Input:
  t - potential Temperature 
  s - Salinity
  optional sigR = xxxx db - reference pressure, default = 2000. db
  allowed t,s as scalars or 1D arrays
  """
# Check type of input var:
  if isinstance(t, float):
    f1d = False
  elif isinstance(t, np.ndarray):
    f1d = True
  else:
    raise Exception(' Type of input t, s should be float or 1d np array')

  if not f1d and s<0.0:
    raise Exception (' Salinity < 0 !!!')
  if f1d and np.min(s)<0.0:
    raise Exception (' Salinity < 0 !!!')

  aone =1.0
  ahalf=1.0/2.0
  a3rd =1.0/3.0
  athird =a3rd
  a4th =1.0/4.0
  afourth=a4th

  """
  c --- Jackett, McDougall, Feistel, Wright and Griffies (2006),
  c --- Algorithms for Density, Potential Temperature, Conservative
  c --- Temperature, and the Freezing Temperature of Seawater, JAOT
  c
  c --- coefficients for 25-term rational function sigloc().
  """
  c001= 9.9984085444849347e+02     #num. constant    coefficent
  c002= 7.3471625860981584e+00     #num.    T        coefficent
  c003=-5.3211231792841769e-02     #num.    T^2      coefficent
  c004= 3.6492439109814549e-04     #num.    T^3      coefficent
  c005= 2.5880571023991390e+00     #num.       S     coefficent
  c006= 6.7168282786692355e-03     #num.    T  S     coefficent
  c007= 1.9203202055760151e-03     #num.       S^2   coefficent
  c008= 1.0000000000000000e+00     #den. constant    coefficent
  c009= 7.2815210113327091e-03     #den.    T        coefficent
  c010=-4.4787265461983921e-05     #den.    T^2      coefficent
  c011= 3.3851002965802430e-07     #den.    T^3      coefficent
  c012= 1.3651202389758572e-10     #den.    T^4      coefficent
  c013= 1.7632126669040377e-03     #den.       S     coefficent
  c014= 8.8066583251206474e-06     #den.    T  S     coefficent
  c015= 1.8832689434804897e-10     #den.    T^3S     coefficent
  c016= 5.7463776745432097e-06     #den.    T  S^1.5 coefficent
  c017= 1.4716275472242334e-09      #den.    T^3S^1.5 coefficent

  c018= 1.1798263740430364e-02     #num. P           coefficent
  c019= 9.8920219266399117e-08     #num. P  T^2      coefficent
  c020= 4.6996642771754730e-06     #num. P     S     coefficent
  c021= 2.5862187075154352e-08     #num. P^2         coefficent
  c022= 3.2921414007960662e-12     #num. P^2T^2      coefficent
  c023= 6.7103246285651894e-06     #den. P           coefficent
  c024= 2.4461698007024582e-17     #den. P^2T^3      coefficent
  c025= 9.1534417604289062e-18      #den. P^3T        coefficent

# --- additional coefficients for dsiglocdt().
  c031= 7.3471625860981580e+00     # num. constant    coefficent
  c032=-1.0642246358568354e-01     # num.    T        coefficent
  c033= 1.0947731732944364e-03     # num.    T^2      coefficent
  c034= 6.7168282786692355e-03     # num.       S     coefficent
  c035= 7.2815210113327090e-03     # den. constant    coefficent
  c036=-8.9574530923967840e-05     # den.    T        coefficent
  c037= 1.0155300889740728e-06     # den.    T^2      coefficent
  c038= 5.4604809559034290e-10     # den.    T^3      coefficent
  c039=-8.8066583251206470e-06     # den.       S     coefficent
  c040= 5.6498068304414700e-10     # den.    T^2S     coefficent
  c041= 2.9432550944484670e-09     # den.    T  S^1.5 coefficent
  c042= 1.9784043853279823e-07     # num. P  T        coefficent
  c043= 6.5842828015921320e-12     # num. P^2T        coefficent
  c044= 7.3385094021073750e-17     # den. P^2T^2      coefficent
  c045= 9.1534417604289060e-18      # den. P^3         coefficent
# --- additional coefficients for dsiglocds().
  c051= 2.5880571023991390e+00     # num. constant    coefficent
  c052= 6.7168282786692355e-03     # num.    T        coefficent
  c053= 3.8406404111520300e-03     # num.       S     coefficent
  c054= 1.7632126669040377e-03     # den. constant    coefficent
  c055=-8.8066583251206470e-06     # den.    T        coefficent
  c056= 1.8832689434804897e-10     # den.    T^3      coefficent
  c057= 8.6195665118148150e-06     # den.       S^0.5 coefficent
  c058= 2.2074413208363504e-09     # den.    T^2S^0.5 coefficent
  c059= 4.6996642771754730e-06     # num. P           coefficent

  sqrmin = 0.0                     # sqrt arg can't be negative
# reference pressure
# sigma0 -> ref pressure = 0 dbar
# sigma2000 -> ref pressure = 2000 dbar
  prs2pdb = 1.e-4                  # Pascals to dbar
  rpdb    = sigR                   # ref. pressure in dbar  
  pref    = rpdb/prs2pdb           # ref. pressure in Pascals, sigma0
# --- coefficients for 17-term rational function sig() at rpdb.  
  c101=c001+(c018-c021*rpdb)*rpdb  # num. constant    coefficent
  c103=c003+(c019-c022*rpdb)*rpdb  # num.    T^2      coefficent
  c105=c005+c020*rpdb              # num.       S     coefficent
  c108=c008+c023*rpdb              # den. constant    coefficent
  c109=c009-c025*rpdb**3           # den.    T        coefficent
  c111=c011-c024*rpdb**2           # den.    T^3      coefficent
# --- additional coefficients for dsigdt().
  c132=c032+(c042-c043*rpdb)*rpdb  # num.    T        coefficent
  c135=c035-c045*rpdb**3           # den. constant    coefficent
  c137=c037-c044*rpdb**2           # den.    T^2      coefficent
# --- additional coefficients for dsigds().
  c151=c051+c059*rpdb              # num. constant    coefficent

  rhoref = 1000.                   # rhoref=qthref kg/m^3

#
# --- -----------------
# --- equation of state
# --- -----------------
#
# --- sigma at rpdb (dbar) as a function of pot.temp (deg c) and salinity (psu)
#
  sqrmin = 0.0
  if f1d:
    s[np.where(s<sqrmin)] = sqrmin

  sig_n = c101 + t*(c002+t*(c103+t*c004)) + \
                         s*(c105-t*c006+s*c007)
  sig_d = c108 + t*(c109+t*(c010+t*(c111+t*c012))) + \
                         s*(c013-t*(c014+t*t*c015) + \
                         np.sqrt(s)*(c016+t*t*c017))
  sig_q = aone/sig_d
  sigma = sig_n*sig_q - rhoref

  return sigma


