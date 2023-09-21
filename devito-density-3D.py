import numpy as np
from devito import (Function, TimeFunction, cos, sin, solve,
                    Eq, Operator, configuration, norm)
from examples.seismic import TimeAxis, RickerSource, Receiver, demo_model
from examples.seismic.model import SeismicModel
from matplotlib import pyplot as plt
from matplotlib import rcParams
import json
rcParams['font.size'] = 12
px = 1/plt.rcParams['figure.dpi']  # pixel in inches

# NBVAL_IGNORE_OUTPUT   
nx = 601
ny = 601
nz = 601
shape   = (nx,ny,nz) 
spacing = (10.,10.,10.) # spacing of 10 meters
origin  = (0.,0.,0.)  
nbl = 0  # number of pad points
sorder=2
torder=2

# # initialize Thomsem parameters to those used in Mu et al., (2020)
vp = 3
rho = 1.0
eps =0.24
delt=0.1

freq_rickers = 0.005
freq = freq_rickers *1# 2.5


model = SeismicModel(space_order=sorder, vp=vp, origin=origin, shape=shape,
                      dtype=np.float32, spacing=spacing, nbl=nbl, epsilon=eps,
                      delta=delt, theta=0., phi=None, bcs="damp")

# Get symbols from model
theta = model.theta
delta = model.delta
epsilon = model.epsilon
m = model.m

# Compute Tmax (just before wave reach border)
vh = vp*np.sqrt(1+2*eps)
vn = vp*np.sqrt(1+2*delt)
# max velocity
vmax = np.maximum(vh,vn)
# dimension length
xmax = ((nx-1) * spacing[0])/2.
xmin = -xmax
ymax = ((ny-1) * spacing[1])/2.
ymin = -ymax
zmax = ((nz-1) * spacing[2])/2.
zmin = -zmax
min_horiz_dist= np.min([np.abs(xmax),np.abs(ymax), np.abs(xmin), np.abs(ymin)])
#Time before wave reaches border
Th = min_horiz_dist / vh # time to reach horizontal border
min_vert_dist= np.minimum(np.abs(zmax),np.abs(zmin))
Tn = min_vert_dist / vn # time to reach top or bottom border
# Choose T max such that wave just reach border (but less than 1 seconds)
Tmax = np.around(np.minimum(Th, Tn), decimals = 3)
Tmax = np.minimum(Tmax, 1000)
# round from dt

# Compute the dt and set time range
t0 = 0.   #  Simulation time start

dt = model.critical_dt
# round Tmax to be a multiple of dt
print('Tmax = ', Tmax)
Tmax = int(Tmax/dt)*(dt) 
Tmax = np.around(Tmax, decimals=4)
ndt_approx = int(Tmax/dt)
print('Tmax = ', Tmax)


# dt = (dvalue/(np.pi*vmax))*np.sqrt(1/(1+etamax*(max_cos_sin)**2)) # eq. above (cell 3)
time_range = TimeAxis(start=t0,stop=Tmax,step=dt)
print("time_range; ", time_range)

# NBVAL_IGNORE_OUTPUT

# time stepping 
p = TimeFunction(name="p", grid=model.grid, time_order=torder, space_order=sorder)#, save=time_range.num)
q = TimeFunction(name="q", grid=model.grid, time_order=torder, space_order=sorder) # space order 4?

# Main equations
# Should use delt instead of delta due to non-float format
term_p = (1 + 2*eps)/rho*(p.dx2 + p.dy2) + np.sqrt(1+2*delt)/rho*q.dz2
term_q = np.sqrt(1 + 2*delt)/rho*(p.dx2 + p.dy2) + (1/rho)*q.dz2

stencil_p = solve(m*p.dt2 - term_p, p.forward)
update_p = Eq(p.forward, stencil_p)

stencil_q = solve(m*q.dt2 - term_q, q.forward)
update_q = Eq(q.forward, stencil_q)

# set source and receivers
src = RickerSource(name='src',grid=model.grid,f0=freq,npoint=1,time_range=time_range)
src.coordinates.data[:,0] = model.domain_size[0]* .5
src.coordinates.data[:,1] = model.domain_size[1]* .5
src.coordinates.data[:,2] = model.domain_size[2]* .5
# Define the source injection
src_term   = src.inject(field=p.forward, expr = -src*dt**2/m )
src_term_q = src.inject(field=q.forward, expr = -src*dt**2/m )

rec  = Receiver(name='rec',grid=model.grid,npoint=shape[0],time_range=time_range)
rec.coordinates.data[:, 0] = np.linspace(model.origin[0],model.domain_size[0], num=model.shape[0])
rec.coordinates.data[:, 1] = 2*spacing[1]
rec.coordinates.data[:, 2] = 2*spacing[1]
# Create interpolation expression for receivers
rec_term = rec.interpolate(expr=p.forward)

# Operators
optime=Operator([update_p] + [update_q] + src_term  + src_term_q+rec_term)

# you can print the generated code for both operators by typing print(optime) and print(oppres)

optime(time=time_range.num-2, dt=model.critical_dt)

# Some useful definitions for plotting if nbl is set to any other value than zero
nxpad,nzpad = shape[0] + 2 * nbl, shape[1] + 2 * nbl
shape_pad   = np.array(shape) + 2 * nbl
origin_pad  = tuple([o - s*nbl for o, s in zip(origin, spacing)])
extent_pad  = tuple([s*(n-1) for s, n in zip(spacing, shape_pad)])


#time_target = 820
ix = int(nx/2)
iy = int(ny/2)
iz = int(nz/2)
ndt = p.data.shape[0]
print("p.data.shape" + str(p.data.shape))
it = 3 #np.minimum(int(time_target/dt), ndt-1)
time = Tmax


###
data = {}
data["nx"] = nx
data["ny"] = ny
data["nz"] = nz
data["xmax"] = xmax
data["xmin"] = xmin
data["ymax"] = ymax
data["ymin"] = ymin
data["zmax"] = zmax
data["zmin"] = zmin

data["ix"] = ix
data["iy"] = iy
data["iz"] = iz


data["dt"] = float(dt)
data["Tmax"] = Tmax
data["epsilon"] = eps
data["delta"] = delt
data["rho"] = rho
data["vp"] = vp
data["freq"] = freq

filename = "devito-density-3D"
print("Saving data on " + filename + ".json")
with open(filename + ".json", 'w', encoding='utf-8') as f:
  json.dump(data, f, ensure_ascii=False, indent=4)

np.save(filename, p.data[2, ::])
print("Saving Numpy arrays on " + filename + ".npy")
