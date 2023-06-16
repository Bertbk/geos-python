import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.size'] = 12
px = 1/plt.rcParams['figure.dpi']  # pixel in inches

# Parser
parser = argparse.ArgumentParser(description='Automatic computation of data (dt, element size, ..) from a JSON file (containing: xmax,xmin,ymax,ymin,zmax,zmin, vp, f, delta, eps, vti_f). Warning: override JSON file!')
parser.add_argument('-f', help='input (root) filename', default="data")
args = parser.parse_args()

# get filename
rootname = args.f

jsonfile = rootname + ".json"
npyfile = rootname + ".npy"

# load data
p = np.load(npyfile)
with open(jsonfile) as f:
  data = json.load(f)

# get data
xmax  = data["xmax"]
xmin  = data["xmin"]
ymax  = data["ymax"]
ymin  = data["ymin"]
zmax  = data["zmax"]
zmin  = data["zmin"]
dt    = data["dt"]
Tmax  = data["Tmax"]
sigma = data["sigma"]
eps   = data["epsilon"]
delt  = data["delta"]
vti_f = data["vti_f"]
vp    = data["vp"]
freq  = data["freq"]

ix = data["ix"]
iy = data["iy"]
iz = data["iz"]
nx = data["nx"]
ny = data["ny"]
nz = data["nz"]

# Plot
time = Tmax

# cut in xz plane
vmax=np.percentile(np.abs(p[:,iy,:]), 99.5)
fig, ax = plt.subplots(figsize=(800*px, 600*px))
pos = ax.imshow(np.transpose(p[:,iy,:]), cmap="seismic", vmin=-vmax, vmax=+vmax, extent=[xmin,xmax, zmin,zmax])
ax.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
ax.set_title("Devito: Wavefield at t="+str(format(time, '.2f'))+"ms with sigma="+str(sigma))
ax.set_xlabel("X Coordinate (m)")
ax.set_ylabel("Z Coordinate (m)")
ax.grid()
fig.colorbar(pos, ax=ax)
figtitle = "devito-fletcher-3D-xz-sigma"+str(sigma) +".png"
print("Saving figure as ... " + figtitle)
fig.savefig(figtitle)
print("done")


# cut in xy plane
vmax=np.percentile(np.abs(p[:,:,iz]), 99.5)
fig, ax = plt.subplots(figsize=(800*px, 600*px))
pos = ax.imshow(np.transpose(p[:,:,iz]), cmap="seismic", vmin=-vmax, vmax=+vmax, extent=[xmin,xmax, ymin,ymax])
ax.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
ax.set_title("Devito: Wavefield at t="+str(format(time, '.2f'))+"ms with sigma="+str(sigma))
ax.set_xlabel("X Coordinate (m)")
ax.set_ylabel("Y Coordinate (m)")
ax.grid()
fig.colorbar(pos, ax=ax)
figtitle = "devito-fletcher-3D-xy-sigma"+str(sigma) +".png"
print("Saving figure as ... " + figtitle)
fig.savefig(figtitle)
print("done")

# cut in yz plane
vmax=np.percentile(np.abs(p[ix,:,:]), 99.5)
fig, ax = plt.subplots(figsize=(800*px, 600*px))
pos = ax.imshow(np.transpose(p[ix,:,:]), cmap="seismic", vmin=-vmax, vmax=+vmax, extent=[ymin,ymax, zmin,zmax])
ax.tick_params('both', length=2, width=0.5, which='major',labelsize=10)
ax.set_title("Devito: Wavefield at t="+str(format(time, '.2f'))+"ms with sigma="+str(sigma))
ax.set_xlabel("Y Coordinate (m)")
ax.set_ylabel("Z Coordinate (m)")
ax.grid()
fig.colorbar(pos, ax=ax)
figtitle = "devito-fletcher-3D-yz-sigma"+str(sigma) +".png"
print("Saving figure as ... " + figtitle)
fig.savefig(figtitle)
print("done")
