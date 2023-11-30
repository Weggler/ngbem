#  testfile for checking hmatrix generation


import sys
sys.path.append('../build')

from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
from libbem import *


#sp = Sphere( (0,0,0), 1)
#mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(1)
mesh = Mesh(unit_cube.GenerateMesh(maxh=1))

#fes = SurfaceL2(mesh, order=1, dual_mapping=True)
fes = H1(mesh, order=3)
u,v = fes.TnT()

a8 = BilinearForm(fes)
SingleLayerPotentialOperator(a8, intorder=8)
#a8.Assemble();
# print (a8.mat)

