#! /usr/bin/python3

# Cut an hdf5 file in three array at the source position
# Save it as json

import argparse
import numpy as np
import json
import h5py
import math
import scipy.ndimage as nd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.widgets import Slider

parser = argparse.ArgumentParser(description='Cut a hdf5 file using JSON data. Warning: override of JSON file')
parser.add_argument('-json', help='input/output json filename', default="data.json")
aarser.add_argument('-i', help='input hdf5 filename', default='pressure_history.hdf5')
args = parser.parse_args()
jsonfile = args.json
hdf5file = args.i

# get data
with open(jsonfile) as f:
  data = json.load(f)

# get results
f=h5py.File(hdf5file)
p=np.array(f['pressure_np1'])
f.close()

print("p old shape = " + str(p.shape))

#get dimension
ndt = p.shape[0] #data["ndt"] might be wrong due to GEOS computing !
nx = data["nx"]
ny = data["ny"]
nz = data["nz"]

p=np.reshape(p,(ndt, nx,ny,nz))
print("p new shape = " + str(p.shape))

# indices to apply the cut (here: on the source point, in the middle)
ix = int(np.floor(nx/2))
iy = int(np.floor(ny/2))
iz = int(np.floor(nz/2))

# Save data to reduce size
np.save("pyz", p[:,ix,:,:])
np.save("pxz", p[:,:,iy,:])
np.save("pxy", p[:,:,:,iz])

# update json data
data["ix"] = ix
data["iy"] = iy
data["iz"] = iz
data["ndt"] = ndt # update of ndt

print("Update " + jsonfile)
with open(jsonfile, 'w', encoding='utf-8') as f:
  json.dump(data, f, ensure_ascii=False, indent=4)