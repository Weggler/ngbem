#  testfile for checking integration
#  commenting out common edge (super fast):
# err-integration (8-10):  2.628425001383119e-06
# err-integration (10-12):  2.402190701995322e-07
# err-integration (12-14):  2.2668739171070717e-08
# err-integration (14-16):  2.1439279138352254e-09
# err-integration (16-18):  2.0177344636389932e-10

#
# all terms (bad):
# err-integration (8-10):  0.0033734366594167387
# err-integration (10-12):  0.001522098631092327
# err-integration (12-14):  0.0007640002683572761
# err-integration (14-16):  0.00041900967880284713
# err-integration (16-18):  0.0002450767472834361

import sys
sys.path.append('../build')

from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
from libbem import *


sp = Sphere( (0,0,0), 1)
mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(1)
# mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))

fes = SurfaceL2(mesh, order=0, dual_mapping=True)
u,v = fes.TnT()

a8 = BilinearForm(fes)
SingleLayerPotential(a8, intorder=8)
a8.Assemble();

a10 = BilinearForm(fes)
SingleLayerPotential(a10, intorder=10)
a10.Assemble();
# print (a.mat)

a12 = BilinearForm(fes)
SingleLayerPotential(a12, intorder=12)
a12.Assemble();

a14 = BilinearForm(fes)
SingleLayerPotential(a14, intorder=14)
a14.Assemble();

a16 = BilinearForm(fes)
SingleLayerPotential(a16, intorder=16)
a16.Assemble();

a18 = BilinearForm(fes)
SingleLayerPotential(a18, intorder=18)
a18.Assemble();


print ("err-integration (8-10): ", Norm (a8.mat.AsVector()-a10.mat.AsVector()))
print ("err-integration (10-12): ", Norm (a10.mat.AsVector()-a12.mat.AsVector()))
print ("err-integration (12-14): ", Norm (a12.mat.AsVector()-a14.mat.AsVector()))
print ("err-integration (14-16): ", Norm (a14.mat.AsVector()-a16.mat.AsVector()))
print ("err-integration (16-18): ", Norm (a16.mat.AsVector()-a18.mat.AsVector()))
