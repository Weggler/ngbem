#  testfile for checking integration

import sys
sys.path.append('../build')

from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
from libbem import *


sp = Sphere( (0,0,0), 1)
mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(1)
# mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))

fes = SurfaceL2(mesh, order=3, dual_mapping=True)
u,v = fes.TnT()

a8 = BilinearForm(fes)
SingleLayerPotentialOperator(a8, intorder=8)
a8.Assemble();

a10 = BilinearForm(fes)
SingleLayerPotentialOperator(a10, intorder=10)
a10.Assemble();
# print (a.mat)

a12 = BilinearForm(fes)
SingleLayerPotentialOperator(a12, intorder=12)
a12.Assemble();

a14 = BilinearForm(fes)
SingleLayerPotentialOperator(a14, intorder=14)
a14.Assemble();

a16 = BilinearForm(fes)
SingleLayerPotentialOperator(a16, intorder=16)
a16.Assemble();

a18 = BilinearForm(fes)
SingleLayerPotentialOperator(a18, intorder=18)
a18.Assemble();


print ("err-integration (8-10): ", Norm (a8.mat.AsVector()-a10.mat.AsVector()))
print ("err-integration (10-12): ", Norm (a10.mat.AsVector()-a12.mat.AsVector()))
print ("err-integration (12-14): ", Norm (a12.mat.AsVector()-a14.mat.AsVector()))
print ("err-integration (14-16): ", Norm (a14.mat.AsVector()-a16.mat.AsVector()))
print ("err-integration (16-18): ", Norm (a16.mat.AsVector()-a18.mat.AsVector()))
