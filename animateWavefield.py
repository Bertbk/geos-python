#! /usr/bin/python3


import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import rcParams
rcParams['font.size'] = 12
px = 1/plt.rcParams['figure.dpi']  # pixel in inches

# get Args
parser = argparse.ArgumentParser(description='Plot, at different time shot, the wavefield of .npy files')
parser.add_argument('-json', help='input json filename', default="data.json")
parser.add_argument('-o', help='output file', default="wavefield")
args = parser.parse_args()
jsonfile = args.json
outputname  = args.o

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
pxz = np.load("pxz.npy")
pyz = np.load("pyz.npy")
pxy = np.load("pxy.npy")
print("pxz.shape = " + str(pxz.shape))
ndt = pxz.shape[0]

#  Time ref (for scale)
time_target = 0.82 # in seconds
itref = np.minimum(int(time_target/dt), ndt-2)
timeref = itref*dt

# Animate function
fps = 30
def animate_func(i):
    im.set_array(pressure[i,:,:])
    ax.set_title("GEOS: Wavefield at t="+str(format(i*dt*1000., '.2f'))+"ms")
    return [im]


# XZ plane
pressure = pxz
vmax=np.percentile(np.abs(pxz[itref,::]), 99.5)
fig, ax = plt.subplots(figsize=(800*px, 600*px))
im = ax.imshow(np.transpose(pxz[1,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[xmin, xmax, zmin, zmax])
ax.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
ax.set_title("GEOS: Wavefield at t="+str(format(dt*1000., '.2f'))+"ms")
ax.set_xlabel("X Coordinate (m)")
ax.set_ylabel("Z Coordinate (m)")
ax.grid()
fig.colorbar(im, ax=ax)
anim = animation.FuncAnimation(
                               fig, 
                               animate_func, 
                               frames = ndt,
                               interval = 200, # in ms
                               )
anim.save('wavefield-xz.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

# XY plane
pressure = pxy
vmax=np.percentile(np.abs(pxy[itref,::]), 99.5)
fig, ax = plt.subplots(figsize=(800*px, 600*px))
im = ax.imshow(np.transpose(pxy[1,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic', extent=[xmin, xmax, ymin, ymax])
ax.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
ax.set_title("GEOS: Wavefield at t="+str(format(dt*1000., '.2f'))+"ms")
ax.set_xlabel("X Coordinate (m)")
ax.set_ylabel("Y Coordinate (m)")
ax.grid()
fig.colorbar(im, ax=ax)
anim = animation.FuncAnimation(
                               fig, 
                               animate_func, 
                               frames = ndt,
                               interval = 200, # in ms
                               )
anim.save('wavefield-xy.mp4', fps=fps, extra_args=['-vcodec', 'libx264'])

print('Done!')


