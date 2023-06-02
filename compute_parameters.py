#! /usr/bin/python3

# exec(open("./makexml.py").read())

import numpy as np
import geos_json

# length (m)
xmax = 3000
xmin = -3000
ymax = 3000
ymin = -3000
zmax = 3000
zmin = -3000

# Source
xs = 0
ys = 0
zs = 0
# frequency of the source (Hz)
f = 20

# pressure velocity (m.s^-1)
vp = 3000

# Thomsen parameters
delta = 0.1
eps = 0.24
vti_f = 1

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

data = geos_json.setData(xmax,xmin,ymax,ymin,zmax,zmin,nx_elem,ny_elem,nz_elem, delta, eps, vti_f,xs,ys,zs, f,cfl_factor, Tmax, dt, ndt, vp,vh,vn, box_eps)
geos_json.write("data.json", data)
