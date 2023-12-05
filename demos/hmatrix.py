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

fesL2 = SurfaceL2(mesh, order=2)
fesH1 = H1(mesh, order=3)
u,v = fesL2.TnT()
uH1, vH1 = fesH1.TnT()

V = BilinearForm(fesL2)
K = BilinearForm(trialspace=fesH1, testspace=fesL2)
SingleLayerPotentialOperator(V, intorder=3)
DoubleLayerPotentialOperator(K, intorder=3)

