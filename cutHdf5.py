#! /usr/bin/python3

# Cut an hdf5 file in three array at the source position
# Save it as json

import numpy as np
import h5py
import json
import math
import scipy.ndimage as nd
import matplotlib.pyplot as plt
from matplotlib import rcParams


# get data

with open('data.json') as f:
  data = json.load(f)

# get results
f=h5py.File('pressure_history.hdf5')
p=np.array(f['pressure_np1'])
f.close()

print(p.shape)

# p=np.reshape(p,(5,188,200,200))
# print(p.shape)

# # matplotlib inline

# dt=0.002
# ix=100
# iy=100
# iz=100
# vmax=0.1

# nt=p.shape[0]
# nx=p.shape[1]
# ny=p.shape[2]
# nz=p.shape[3]
# ixmin=15
# ixmax=nx-1-ixmin
# iymin=15
# iymax=ny-1-iymin
# izmin=15
# izmax=nz-1-izmin

# def plotWfld(d,it,output="figure",save=False):

#     time=it*dt
#     plt.figure(figsize=(16,6))
#     plt.subplot(1,3,1)
#     plt.imshow(np.transpose(p[it,ix,:,:]),vmin=-vmax,vmax=vmax,cmap='seismic')
#     #plt.hlines((izmin,izmax),iymin,iymax,LineStyle='--')
#     #plt.vlines((iymin,iymax),izmin,izmax,LineStyle='--')
#     plt.xlabel("Y")
#     plt.ylabel("Z")
#     plt.title(r"X-slice @ %0.2f sec" %time)

#     plt.subplot(1,3,2)
#     plt.imshow(np.transpose(p[it,:,iy,:]),vmin=-vmax,vmax=vmax,cmap='seismic')
#     #plt.hlines((izmin,izmax),ixmin,ixmax,LineStyle='--')
#     #plt.vlines((ixmin,ixmax),izmin,izmax,LineStyle='--')
#     plt.xlabel("X")
#     plt.ylabel("Z")
#     plt.title(r"Y-slice @ %0.2f sec" %time)

#     plt.subplot(1,3,3)
#     plt.imshow(np.transpose(p[it,:,:,iz]),vmin=-vmax,vmax=vmax,cmap='seismic')
#     #plt.hlines((iymin,iymax),ixmin,ixmax,LineStyle='--')
#     #plt.vlines((ixmin,ixmax),iymin,iymax,LineStyle='--')
#     plt.xlabel("X")
#     plt.ylabel("Y")
#     plt.title(r"Z-slice @ %0.2f sec" %time)

#     plt.tight_layout()
#     if save==True:
#         plt.savefig("./fig/"+output+"_%d" %it)