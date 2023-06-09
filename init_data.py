#! /usr/bin/python3

import argparse
import numpy as np
import json

parser = argparse.ArgumentParser(description='Init a JSON file')
parser.add_argument('-o', help='output json filename', default="data.json")
args = parser.parse_args()
output_file = args.o

data = {}
# length (m)
data["xmax"] = 3000
data["xmin"] = -3000
data["ymax"] = 3000
data["ymin"] = -3000
data["zmax"] = 3000
data["zmin"] = -3000

# Source
data["xs"] = 0
data["ys"] = 0
data["zs"] = 0
# frequency of the source (Hz)
data["f"] = 20

# pressure velocity (m.s^-1)
data["vp"] = 3000

# Thomsen parameters
data["delta"] = 0.1
data["eps"] = 0.24
data["sigma"] = -1
data["vs"]    = 0
data["vti_f"] = 1

data["nlambda"] = 15 # number of point per wavelength

data["nx_elem"] = -1
data["ny_elem"] = -1
data["nz_elem"] = -1
data["cfl_factor"] = -1
data["Tmax"] = -1
data["dt"] = -1
data["vh"] = -1
data["vn"] = -1
data["box_eps"] =-1
data["dt_hdf5"] = -1
data["dt_vtk"] = -1


with open(output_file, 'w', encoding='utf-8') as f:
  json.dump(data, f, ensure_ascii=False, indent=4)
