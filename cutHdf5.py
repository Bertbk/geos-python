#! /usr/bin/python3

# Cut an hdf5 file in three array at the source position
# Save it as json

import numpy as np
import h5py
import geos_json
import math
import scipy.ndimage as nd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.widgets import Slider

# get data
data = geos_json.getData("data.json")

# get results
f=h5py.File('pressure_history.hdf5')
p=np.array(f['pressure_np1'])
f.close()

print(p.shape)

ndt = 1000#data["ndt"]
nx = 51#data["nx_elem"] +1
ny = 51#data["ny_elem"] +1
nz = 51#data["nz_elem"] +1

p=np.reshape(p,(ndt, nx,ny,nz))
print(p.shape)

# indices to apply the cut (here: on the source point, in the middle)
ix = int(np.floor(nx/2))
iy = int(np.floor(ny/2))
iz = int(np.floor(nz/2))

# Save data to reduce size

np.save("pyz", p[:,ix,:,:])
np.save("pxz", p[:,:,iy,:])
np.save("pxy", p[:,:,:,iz])

# pyz = p[:,ix,:,:]
# pxz = p[:,:,iy,:]
# pxy = p[:,:,:,iz]

# update json data

data["ix"] = ix
data["iy"] = iy
data["iz"] = iz
data["nx"] = nx
data["ny"] = ny
data["nz"] = nz
data["ndt"] = ndt

geos_json.write("data.json", data)