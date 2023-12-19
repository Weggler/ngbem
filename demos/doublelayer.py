#  testfile for checking integration 
#  output of testrun 29.11.23
# err-integration (8-10):  1.6205859881425085e-07
# err-integration (10-12):  1.5225443842172683e-08
# err-integration (12-14):  1.3555959613067985e-09
# err-integration (14-16):  1.3309124393862712e-10
# err-integration (16-18):  1.507548088767588e-11

import sys
sys.path.append('../build')

from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
from libbem import *


sp = Sphere( (0,0,0), 1)
mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(1)
# mesh = Mesh(unit_cube.GenerateMesh(maxh=0.3))

fes = SurfaceL2(mesh, order=0, dual_mapping=False)
fesh1 = H1(mesh, order=1)
u,v = fes.TnT()
uH1, vH1 = fesh1.TnT()

b8 = BilinearForm(trialspace=fesh1, testspace=fes)
DoubleLayerPotentialOperator(b8, intorder=8)
b8.Assemble();

b10 = BilinearForm(trialspace=fesh1, testspace=fes)
DoubleLayerPotentialOperator(b10, intorder=10)
b10.Assemble();
# print (b.mat)

b12 = BilinearForm(trialspace=fesh1, testspace=fes)
DoubleLayerPotentialOperator(b12, intorder=12)
b12.Assemble();

b14 = BilinearForm(trialspace=fesh1, testspace=fes)
DoubleLayerPotentialOperator(b14, intorder=14)
b14.Assemble();

b16 = BilinearForm(trialspace=fesh1, testspace=fes)
DoubleLayerPotentialOperator(b16, intorder=16)
b16.Assemble()

b18 = BilinearForm(trialspace=fesh1, testspace=fes)
DoubleLayerPotentialOperator(b18, intorder=18)
b18.Assemble();


print ("err-integration (8-10): ", Norm (b8.mat.AsVector()-b10.mat.AsVector()))
print ("err-integration (10-12): ", Norm (b10.mat.AsVector()-b12.mat.AsVector()))
print ("err-integration (12-14): ", Norm (b12.mat.AsVector()-b14.mat.AsVector()))
print ("err-integration (14-16): ", Norm (b14.mat.AsVector()-b16.mat.AsVector()))
print ("err-integration (16-18): ", Norm (b16.mat.AsVector()-b18.mat.AsVector()))
