from igakit.io import PetIGA,VTK
from numpy import linspace
import glob
import os

nrb = PetIGA().read("igasol.dat")
#uniform = lambda U: linspace(U[0], U[-1], 400)
for infile in glob.glob("sol*.dat"):
	name = infile.split(".")[0]
	number = name.split("l")[1]
	root = './solV'+number+'.vtk'
	if not os.path.isfile(root):
		sol = PetIGA().read_vec(infile,nrb)
		outfile = './vtkOut/solV'+number+'.vtk' 
		VTK().write(outfile,  
	            nrb,             
	            fields=sol,     
#	            sampler=uniform, 
    		    scalars = {'ice':0,'air':1,'temp':2,'rhov':3}) #,'prDarc':2},
#		    gradients ={'presgrad':[2,2]})
#	    	    vectors = {'vel_darc':[3,4]})
