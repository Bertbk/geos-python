#! /usr/bin/python3

# Cut an hdf5 file in three array at the source position
# Save it as json

import numpy as np
import geos_json
import scipy.ndimage as nd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.widgets import Slider

# get data
data = geos_json.getData("data.json")

ndt = data["ndt"]
dt = data["dt"]
nx = data["nx"]
ny = data["ny"]
nz = data["nz"]

# load pressure

pyz = np.load("pyz.npy")
pxz = np.load("pxz.npy")
pxy = np.load("pxy.npy")

# Matplotlib inline

# indices to apply the cut (here: on the source point, in the middle)
ix = int(np.floor(nx/2))
iy = int(np.floor(ny/2))
iz = int(np.floor(nz/2))
# max value (for scaling)
value_max = np.percentile(np.abs(pxy),99.2)

ixmin = 1
ixmax = nx-1
iymin = 1
iymax = ny-1
izmin = 1
izmax = nz-1

# Function for widget
def plotWfld(dyz,dxz,dxy,fig, ax, it,output="figure",save=False):
  time=it*dt
  ax[0].imshow(np.transpose(dyz[it,:,:]),vmin=-value_max,vmax=value_max,cmap='seismic')
  ax[0].grid()
  ax[0].set_xlabel("Y")
  ax[0].set_ylabel("Z")
  ax[0].set_title(r"X-slice @ %0.2f sec" %time)
  ax[1].imshow(np.transpose(dxz[it,:,:]),vmin=-value_max,vmax=value_max,cmap='seismic')
  ax[1].grid()
  ax[1].set_xlabel("X")
  ax[1].set_ylabel("Z")
  ax[1].set_title(r"Y-slice @ %0.2f sec" %time)
  ax[2].imshow(np.transpose(dxy[it,:,:]),vmin=-value_max,vmax=value_max,cmap='seismic')
  ax[2].grid()
  ax[2].set_xlabel("X")
  ax[2].set_ylabel("Y")
  ax[2].set_title(r"Z-slice @ %0.2f sec" %time)

  fig.tight_layout()
  if save==True:
      plt.savefig("./fig/"+output+"_%d" %it)

#plt.figure(figsize=(16,6))
fig, ax = plt.subplots(1,3)
fig.set_figheight(8)
fig.set_figwidth(16)


# Slider
ax_time = fig.add_axes([0.25, 0.1, 0.65, 0.03])

time_slide = Slider(
  ax_time, "Iteration", 0, ndt,
  valinit=1, valstep=1,
  initcolor='none'  # Remove the line marking the valinit position.
)


def update(val):
    current_time = time_slide.val
    plotWfld(pyz,pxz,pxy,fig,ax,current_time,output="figure",save=False)
    fig.canvas.draw_idle()

plotWfld(pyz,pxz,pxy,fig, ax, int(ndt/2),output="figure",save=False)
time_slide.on_changed(update)
plt.show()
