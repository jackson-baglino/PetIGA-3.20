import numpy as np
import glob
import os
from igakit.io import PetIGA, VTK
from numpy import linspace

nrb = PetIGA().read("igasol.dat")
uniform = lambda U: linspace(U[0], U[-1], 480)

# Import ice grains:
for infile in glob.glob("sol*.dat"):
    name = infile.split(".")[0]
    number = name.split("l")[1]
    root = './solV'+number+'.vtk'

    sol = PetIGA().read_vec(infile,nrb)
    if os.path.exists(root):
        os.remove(root)
    VTK().write(root, 
                nrb, fields=sol, 
                sampler=uniform, 
                scalars={'IcePhase':0, 'Temperature':1, 'VaporDensity':2})
    print("Created: " + root + ".\n")

# Import sediment grains:
nrb2 = None
if os.path.exists("igasoil.dat"):
    nrb2 = PetIGA().read("igasoil.dat")

    for infile in glob.glob("soil*.dat"):
        name = infile.split(".")[0]
        number = name.split("l")[1]
        root = './soilV'+number+'.vtk'

        sol = PetIGA().read_vec(infile,nrb2)
        if os.path.exists(root):
            os.remove(root)
        VTK().write(root, 
                    nrb2, fields=sol, 
                    sampler=uniform, 
                    scalars={'SedPhase':0})
        print("Created: " + root + ".\n")

