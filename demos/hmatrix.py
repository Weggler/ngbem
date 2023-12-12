#  testfile for checking hmatrix generation


import sys
sys.path.append('../build')

from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
from libbem import *


#sp = Sphere( (0,0,0), 1)
#mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(1)
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.1))

fesL2 = SurfaceL2(mesh, order=0)
fesH1 = H1(mesh, order=1, definedon=mesh.Boundaries(".*"))
u,v = fesL2.TnT()
uH1, vH1 = fesH1.TnT()

V = BilinearForm(fesL2)
#K = BilinearForm(trialspace=fesH1, testspace=fesL2)
SingleLayerPotentialOperator(V, intorder=5, leafsize=40, eta=3., eps=1e-4, method="svd")
#DoubleLayerPotentialOperator(K, intorder=3)

