#! /usr/bin/python3

# Cut an hdf5 file in three array at the source position
# Save it as json


import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.size'] = 12
px = 1/plt.rcParams['figure.dpi']  # pixel in inches

# get Args
parser = argparse.ArgumentParser(description='Cut a hdf5 file using JSON data. Warning: override of JSON file')
parser.add_argument('-json', help='input json filename', default="data.json")
args = parser.parse_args()
jsonfile = args.json

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

sigma = data["sigma"]

nx = data["nx_elem"] + 1
ny = data["ny_elem"] + 1
nz = data["nz_elem"] + 1

dt=data["dt_hdf5"]

xscale = np.linspace(xmin, xmax, num =nx, endpoint=True)
yscale = np.linspace(ymin, ymax, num =ny, endpoint=True)
zscale = np.linspace(zmin, zmax, num =nz, endpoint=True)

# load pressure files
pyz = np.load("pyz.npy")
pxz = np.load("pxz.npy")
pxy = np.load("pxy.npy")
print("pyz.shape = " + str(pyz.shape))
ndt = pyz.shape[0]

# Save fig
time_target = 0.82 # in seconds
it = np.minimum(int(time_target/dt), ndt-2)
time = it*dt

# XZ plane
vmax=np.percentile(np.abs(pxz[it,::]), 99.5)
fig, ax = plt.subplots(figsize=(800*px, 600*px))
pos = ax.imshow(np.transpose(pxz[it,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[xmin, xmax, zmin, zmax])
ax.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
ax.set_title("GEOS: Wavefield at t="+str(format(time*1000., '.2f'))+"ms with sigma="+str(sigma))
ax.set_xlabel("X Coordinate (m)")
ax.set_ylabel("Z Coordinate (m)")
ax.grid()
fig.colorbar(pos, ax=ax)

figtitle = "geos-sigma-xz-"+str(sigma)+".png"
print("Saving figure as ... " + figtitle)
fig.savefig(figtitle)
print("done xz" )

# YZ plane
vmax=np.percentile(np.abs(pyz[it,::]), 99.5)
fig, ax = plt.subplots(figsize=(800*px, 600*px))
pos = ax.imshow(np.transpose(pyz[it,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[ymin, ymax, zmin, zmax])
ax.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
ax.set_title("GEOS: Wavefield at t="+str(format(time*1000., '.2f'))+"ms with sigma="+str(sigma))
ax.set_xlabel("Y Coordinate (m)")
ax.set_ylabel("Z Coordinate (m)")
ax.grid()
fig.colorbar(pos, ax=ax)

figtitle = "geos-sigma-yz-"+str(sigma)+".png"
print("Saving figure as ... " + figtitle)
fig.savefig(figtitle)
print("done xy")

# XY plane
vmax=np.percentile(np.abs(pxy[it,::]), 99.5)
fig, ax = plt.subplots(figsize=(800*px, 600*px))
pos = ax.imshow(np.transpose(pxy[it,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[xmin, xmax, ymin, ymax])
ax.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
ax.set_title("GEOS: Wavefield at t="+str(format(time*1000., '.2f'))+"ms with sigma="+str(sigma))
ax.set_xlabel("X Coordinate (m)")
ax.set_ylabel("Y Coordinate (m)")
ax.grid()
fig.colorbar(pos, ax=ax)

figtitle = "geos-sigma-xy-"+str(sigma)+".png"
print("Saving figure as ... " + figtitle)
fig.savefig(figtitle)
print("done xy")