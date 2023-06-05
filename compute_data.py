#! /usr/bin/python3

import argparse
import numpy as np
import json

parser = argparse.ArgumentParser(description='Automatic computation of data (dt, element size, ..) from a JSON file (containing: xmax,xmin,ymax,ymin,zmax,zmin, vp, f, delta, eps, vti_f). Warning: override JSON file!')
parser.add_argument('-o', help='input and output json filename', default="data.json", required=True)
args = parser.parse_args()
filename = args.o

with open(filename) as f:
  data = json.load(f)

# length (m)
xmax = data["xmax"]
xmin = data["xmin"]
ymax = data["ymax"]
ymin = data["ymin"]
zmax = data["zmax"]
zmin = data["zmin"]
# frequency of the source (Hz)
f = data["f"]
# pressure velocity (m.s^-1)
vp = data["vp"]
# Thomsen parameters
delta = data["delta"]
eps   = data["eps"]
vti_f = data["vti_f"]

# Compute size of element
# Horizontal vh and vertical vn velocity
vh = vp*np.sqrt(1+2*eps)
vn = vp*np.sqrt(1+2*delta)

# max velocity
vmax = np.maximum(vh,vn)

# Tmax
# Automatic (when wave touch borders)
min_horiz_dist= np.min([np.abs(xmax),np.abs(ymax), np.abs(xmin), np.abs(ymin)])
Th = min_horiz_dist / vh # time to reach horizontal border
min_vert_dist= np.minimum(np.abs(zmax),np.abs(zmin))
Tn = min_vert_dist / vn # time to reach top or bottom border

# Choose T max such that wave just reach border (but less than 1 seconds)
Tmax = np.around(np.minimum(Th, Tn), decimals = 3)
Tmax = np.minimum(Tmax, 1)
print('Tmax = ', Tmax)

# Wavenumber and wavelength
omega = 2*np.pi*f
k = omega / vmax
wavelength = 2*np.pi/k

# CFL, Space and Time steps (Dx and Dt)
cfl_factor = 0.25
dx = wavelength / 10
dt = np.around(cfl_factor*dx/vmax, decimals = 4)
ndt = int(Tmax / dt)

# Number of hexa in each dimension
nx_elem = int((xmax-xmin)/dx)
ny_elem = int((ymax-ymin)/dx)
nz_elem = int((zmax-zmin)/dx)

box_eps = np.minimum(0.1, np.around(dx/10, decimals = 1))

print("Summary :")
print("Tmax = " + str(Tmax))
print("dt = " + str(dt))
print("ndt = " + str(ndt))
print("nx_elem = " + str(nx_elem))
print("ny_elem = " + str(ny_elem))
print("nz_elem = " + str(nz_elem))

# update data
data["cfl_factor"] = cfl_factor 
data["dx"]         = dx
data["dt"]         = dt 
data["ndt"]        = ndt 
data["nx_elem"]    = nx_elem 
data["ny_elem"]    = ny_elem 
data["nz_elem"]    = nz_elem
data["nx"]         = nx_elem +1
data["ny"]         = ny_elem +1
data["nz"]         = nz_elem +1
data["box_eps"]    = box_eps 
data["omega"]      = omega
data["wavenumber"] = k
data["wavelength"] = wavelength
data["Tmax"]       = Tmax

# write data

with open(output_file, 'w', encoding='utf-8') as f:
  json.dump(data, f, ensure_ascii=False, indent=4)
