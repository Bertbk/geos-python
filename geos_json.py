import json

def setData(xmax,xmin,ymax,ymin,zmax,zmin,nx_elem,ny_elem,nz_elem, delta, eps, vti_f,xs,ys,zs, f,cfl_factor, Tmax, dt, ndt, vp,vh,vn,box_eps):
  data = {}
  data["xmax"] = xmax
  data["xmin"] = xmin
  data["ymax"] = ymax
  data["ymin"] = ymin
  data["zmax"] = zmax
  data["zmin"] = zmin
  data["nx_elem"] = nx_elem
  data["ny_elem"] = ny_elem
  data["nz_elem"] = nz_elem
  data["delta"] = delta
  data["eps"] = eps
  data["vti_f"] = vti_f
  data["xs"] = xs
  data["ys"] = ys
  data["zs"] = zs
  data["f"] = f
  data["cfl_factor"] = cfl_factor
  data["Tmax"] = Tmax
  data["dt"] = dt
  data["ndt"] = ndt
  data["vp"] = vp
  data["vh"] = vh
  data["vn"] = vn
  data["box_eps"] =box_eps
  return data

def write(filename, data):
  with open(filename, 'w', encoding='utf-8') as f:
      json.dump(data, f, ensure_ascii=False, indent=4)

def getData(filename):
  with open(filename) as f:
    data = json.load(f)
  return data