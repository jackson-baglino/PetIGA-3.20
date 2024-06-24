import numpy as np
import glob
import os
from igakit.io import PetIGA, VTK
from numpy import linspace

Nx = os.getenv("Nx")
Ny = os.getenv("Ny")
Nz = os.getenv("Nz")

if Nx is not None and Ny is not None and Nz is not None:
    num_points = max(int(Nx), int(Ny), int(Nz))
else:
    print("Nx, Ny, Nz not set, using default value of 400.")
    num_points = 400

nrb = PetIGA().read("igasol.dat")
uniform = lambda U: linspace(U[0], U[-1], num_points)

# Import ice grains:
counter = 0

for infile in glob.glob("sol*.dat"):
    name = infile.split(".")[0]
    number = name.split("l")[1]
    root = './vtkOut/solV_'+str(number)+'.vtk'

    sol = PetIGA().read_vec(infile,nrb)
    if os.path.exists(root):
        os.remove(root)
    VTK().write(root,
                nrb, fields=sol, 
                sampler=uniform, 
                scalars={'IcePhase':0, 'Temperature':1, 'VaporDensity':2})
    print("Created: " + root + ".\n")

    counter = counter + 1

# Import sediment grains:
nrb2 = None
if os.path.exists("igasoil.dat"):
    nrb2 = PetIGA().read("igasoil.dat")

    for infile in glob.glob("soil*.dat"):
        name = infile.split(".")[0]
        number = name.split("l")[1]
        root = './vtkOut/soilV'+number+'.vtk'

        sol = PetIGA().read_vec(infile,nrb2)
        if os.path.exists(root):
            os.remove(root)
        VTK().write(root, 
                    nrb2, fields=sol, 
                    sampler=uniform, 
                    scalars={'SedPhase':0})
        print("Created: " + root + ".\n")

