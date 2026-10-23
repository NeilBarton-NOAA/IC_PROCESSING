"""
  Functions for plotting and figure settings
"""
import os
import numpy as np
import matplotlib.pyplot as plt

def bottom_text(btx, gtdir='RTOFS_utilities', \
                ipwd=1, pos=[0.05, 0.05], fsz=8, f_short=False):
  sgit = 'github.com/DmitryDukhovskoy/' + gtdir
  drr = os.getcwd()
#  btnm = drr+'/'+btx
  bnm = os.path.basename(drr)
  if f_short:
    btnm = bnm + '/' + btx
  else:
    btnm = sgit + '/' + bnm + '/' + btx

  if ipwd==0:
    btnm = btx

  plt.text(pos[0],pos[1],btnm,horizontalalignment='left',
         transform=plt.gcf().transFigure, fontsize=fsz)
  return


def correct_evensp_grid(x, y):
  xcrct = x.copy()
  ycrct = y.copy()
  if x.ndim == 1:
    pass
  elif x.ndim == 2:
    x_row = x[0, :]
    if not np.allclose(x_row, x):
      raise ValueError("The rows of 'x' must be equal")
  else:
    raise ValueError("'x' can have at maximum 2 dimensions")
#   
  if y.ndim == 1:
    pass
  elif y.ndim == 2:
    y_col = y[:, 0]
    if not np.allclose(y_col, y.T):
      raise ValueError("The columns of 'y' must be equal")
  else:
    raise ValueError("'y' can have at maximum 2 dimensions")

  nx = len(x_row)
  ny = len(y_col)

  width = x_row[-1] - x_row[0]
  height = y_col[-1] - y_col[0]

  rtol = 1.e-6
  atol = 0.0
  dx = width/(nx-1)
  if not np.allclose(np.diff(x_row), dx, atol, rtol):
    x1r = np.arange(x_row[0],x_row[-1]+0.1*dx,dx)

    for ii in range(ny):
      xcrct[ii] = x1r

# breakpoint() 
  dy = height/(ny-1)
  if not np.allclose(np.diff(y_col), dy, atol, rtol):
    y1r = np.arange(y_col[0],y_col[-1]+0.1*dy,dy)

    for ii in range(nx):
      ycrct[:,ii] = y1r

# breakpoint()
  x_row2 = xcrct[0, :]
  np.allclose(np.diff(x_row2), dx, atol, rtol)

  return xcrct, ycrct

def plot_fld2D(aa,ctitle,cl1,cl2,X=np.array([]),Y=np.array([]),clrmp='turbo',
               btx='mod_utils_fig.py', figout='plot_field.jpg'):
  """
    Create a quick figure of 2D map and 
    put it on web to view via browser
    otherwise python get slow
    open figure in the browser:
        https://www-dev.star1.nesdis.noaa.gov/ncei_wod/pthn_fig.jpg
  """
  plt.ion()
  cmpS = plt.cm.get_cmap(clrmp)
  fig1 = plt.figure(1,figsize=(8,8), constrained_layout=False)
  plt.clf()
  ax1 = plt.axes([0.1, 0.1, 0.8, 0.8])
  jdm=aa.shape[0]
  idm=aa.shape[1]

#  if np.isfinite(X) and np.isfinite(Y):
  if X.size == 0 | Y.size == 0:
    im = ax1.pcolormesh(X,Y,aa,cmap=cmpS,shading='flat')
    ax1.axis('scaled')
    ax1.set(xlim=(np.min(X),np.max(X)),
            ylim=(np.min(Y),np.max(Y)))
  else:
    im = ax1.pcolormesh(aa,cmap=cmpS,shading='flat')
    ax1.axis('scaled')
    ax1.set(xlim=(0,idm),ylim=(0,jdm))

  im.set_clim(cl1,cl2)
  ax1.autoscale(enable=True, axis='both', tight=True)

  ax1.set_title(ctitle)
#
# Colorbar
  divider = make_axes_locatable(ax1)
  cax = divider.append_axes("bottom",size="5%",pad=0.25)
  fig1.colorbar(im,ax=ax1,orientation='horizontal',cax=cax)


  bottom_text(btx)

  plt.show()
  plt.savefig(figout)
  plt.close(fig1)

  wwwpth='/net/www-dev/www/ncei_wod/'
  wwwfig='pthn_fig.jpg'
  cmd = ('mv {0} {1}{2}'.format(figout,wwwpth,wwwfig))
  os.system(cmd)


  return


