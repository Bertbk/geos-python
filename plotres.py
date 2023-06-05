#! /usr/bin/python3

# Cut an hdf5 file in three array at the source position
# Save it as json


import argparse
import json
import numpy as np
import scipy.ndimage as nd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.widgets import Slider

# get data
parser = argparse.ArgumentParser(description='Cut a hdf5 file using JSON data. Warning: override of JSON file')
parser.add_argument('-json', help='input json filename', default="data.json")
args = parser.parse_args()
jsonfile = args.json

# get data
with open(jsonfile) as f:
  data = json.load(f)

ndt = data["ndt"]
dt = data["dt"]
nx = data["nx"]
ny = data["ny"]
nz = data["nz"]

xmin = data["xmin"]
xmax = data["xmax"]
ymin = data["ymin"]
ymax = data["ymax"]
zmin = data["zmin"]
zmax = data["zmax"]

xscale = np.linspace(xmin, xmax, num =nx, endpoint=True)
yscale = np.linspace(ymin, ymax, num =ny, endpoint=True)
zscale = np.linspace(zmin, zmax, num =nz, endpoint=True)

# load pressure files
pyz = np.load("pyz.npy")
pxz = np.load("pxz.npy")
pxy = np.load("pxy.npy")

# Matplotlib inline

# max value (for scaling)
value_max = np.percentile(np.abs(pxy),99.2)

# ixmin = 1
# ixmax = nx-1
# iymin = 1
# iymax = ny-1
# izmin = 1
# izmax = nz-1

# Function for widget
def plotWfld(dyz,dxz,dxy,fig, ax, it,output="figure",save=False):
  time=it*dt
  ax[0].imshow(np.transpose(dyz[it,:,:]),vmin=-value_max,vmax=value_max,cmap='seismic')
  ax[0].grid()
  ax[0].set_xlabel("Y")
  ax[0].set_ylabel("Z")
  ax[0].set_title(r"X-slice @ %0.2f sec" %time)
  ax[0].set_xlim([ymin, ymax])
  ax[0].set_ylim([zmin, zmax])
  ax[1].imshow(np.transpose(dxz[it,:,:]),vmin=-value_max,vmax=value_max,cmap='seismic')
  ax[1].grid()
  ax[1].set_xlabel("X")
  ax[1].set_ylabel("Z")
  ax[1].set_title(r"Y-slice @ %0.2f sec" %time)
  ax[1].set_xlim([xmin, xmax])
  ax[1].set_ylim([zmin, zmax])
  ax[2].imshow(np.transpose(dxy[it,:,:]),vmin=-value_max,vmax=value_max,cmap='seismic')
  ax[2].grid()
  ax[2].set_xlabel("X")
  ax[2].set_ylabel("Y")
  ax[2].set_title(r"Z-slice @ %0.2f sec" %time)
  ax[2].set_xlim([xmin, xmax])
  ax[2].set_ylim([ymin, ymax])

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
