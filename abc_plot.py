#! /usr/bin/python3

# Cut an hdf5 file in three array at the source position
# Save it as json


import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.animation as animation


rcParams['font.size'] = 12
px = 1/plt.rcParams['figure.dpi']  # pixel in inches

# get Args
parser = argparse.ArgumentParser(description='Plot results for VTI assuming .npy files exist')
parser.add_argument('-json', help='input json filename', default="data.json")
parser.add_argument('-t', help='Start Time', default="0.8")
parser.add_argument('-n', help='Number of screenshot', default="1")
parser.add_argument('-animated', help='Animation? (boolean)', default="False")
args = parser.parse_args()
jsonfile = args.json
time_start = float(args.t)
is_animated = bool(args.animated)
nshots = int(args.n)

# get data
with open(jsonfile, mode="r") as f:
  data = json.load(f)

dt = data["dt_hdf5"]
nx = data["nx_elem"] + 1
ny = data["ny_elem"] + 1
nz = data["nz_elem"] + 1

xmin = data["xmin"]
xmax = data["xmax"]
ymin = data["ymin"]
ymax = data["ymax"]
zmin = data["zmin"]
zmax = data["zmax"]

nx = data["nx_elem"] + 1
ny = data["ny_elem"] + 1
nz = data["nz_elem"] + 1

dt=data["dt_hdf5"]

xscale = np.linspace(xmin, xmax, num =nx, endpoint=True)
yscale = np.linspace(ymin, ymax, num =ny, endpoint=True)
zscale = np.linspace(zmin, zmax, num =nz, endpoint=True)

# load pressure files
#pyz = np.load("pyz.npy")
pxz_vti = np.load("pxz_vti.npy")
pxy_vti = np.load("pxy_vti.npy")

pxz_iso = np.load("pxz_iso.npy")
pxy_iso = np.load("pxy_iso.npy")

print("pxz.shape = " + str(pxz_vti.shape))
ndt = pxz_vti.shape[0]

# Save fig
it_start = np.minimum(int(time_start/dt), ndt-2)
time_start = it_start*dt
if(it_start == ndt-2):
  print("Time start is too large, time = "+time_start)
else:
  print("Time start= " + str(time_start))

if( nshots > (ndt - it_start)):
  nshots = ndt - it_start

tshots = np.linspace(it_start, ndt - 2, 10, dtype=int)

if(is_animated):
  # lists of lists (of artistic)
  ims = []


fig, ((axIso_xz, axVti_xz), (axIso_xy, axVti_xy))  = plt.subplots(2,2, figsize=(1600*px, 1200*px))
axIso_xz.set_xlabel("X Coordinate (m)")
axIso_xz.set_ylabel("Z Coordinate (m)")
axIso_xz.grid()
axIso_xy.set_xlabel("X Coordinate (m)")
axIso_xy.set_ylabel("Y Coordinate (m)")
axIso_xy.grid()
axVti_xz.set_xlabel("X Coordinate (m)")
axVti_xz.set_ylabel("Z Coordinate (m)")
axVti_xz.grid()
axVti_xy.set_xlabel("X Coordinate (m)")
axVti_xy.set_ylabel("Y Coordinate (m)")
axVti_xy.grid()

axIso_xz.set_title("Isotropic ABC")
axVti_xz.set_title("Tuned ABC")
# XZ plane
vmaxIso=np.percentile(np.abs(pxz_iso[tshots[0],::]), 99.5)
vmaxVti=np.percentile(np.abs(pxz_vti[tshots[0],::]), 99.5)
vmax = vmaxIso

for idit, it in enumerate(tshots):
  time = it*dt
  print("Printing fig for time=" +str(format(time*1000., '.2f'))+"ms")
  fig.suptitle("GEOS: Wavefield at t="+str(format(time*1000., '.2f'))+"ms")

  axIso_xz.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
  posIso_xz = axIso_xz.imshow(np.transpose(pxz_iso[it,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[xmin, xmax, zmin, zmax], animated=True )
  if(idit == 0):
    fig.colorbar(posIso_xz, ax=axIso_xz)

  axVti_xz.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
  posVti_xz = axVti_xz.imshow(np.transpose(pxz_vti[it,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[xmin, xmax, zmin, zmax], animated=True )
  if(idit == 0):
    fig.colorbar(posVti_xz, ax=axVti_xz)

  axIso_xy.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
  posIso_xy = axIso_xy.imshow(np.transpose(pxy_iso[it,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[xmin, xmax, zmin, zmax], animated=True )
  if(idit == 0):
    fig.colorbar(posIso_xy, ax=axIso_xy)

  axVti_xy.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
  posVti_xy = axVti_xy.imshow(np.transpose(pxy_vti[it,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[xmin, xmax, zmin, zmax], animated=True )
  if(idit == 0):
    fig.colorbar(posVti_xy, ax=axVti_xy)

  figtitle = "ABC-comparison-t-"+str(np.floor(time*1000))+".png"
  print("Saving figure as ... " + figtitle)
  fig.savefig(figtitle)

  if(is_animated):
    ims.append([posIso_xz, posVti_xz, posIso_xy, posVti_xy])

# ========================== Animate
# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame


if(is_animated):
  ani = animation.ArtistAnimation(fig, ims, interval=500, blit=True, repeat_delay=1000, repeat = True)
  ani.save('mymovie.mp4')
  
