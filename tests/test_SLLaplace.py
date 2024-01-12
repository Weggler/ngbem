from netgen.occ import *
from ngsolve import *
from ngsolve.webgui import Draw
from libbem import *
from ngsolve import Projector, Preconditioner
from ngsolve.krylovspace import CG


sp = Sphere( (0,0,0), 1)
mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.3)).Curve(3)

fesL2 = SurfaceL2(mesh, order=3, dual_mapping=True)
u,v = fesL2.TnT()

V = SingleLayerPotentialOperator(fesL2, intorder=10, method="svd", leafsize=40, eta=3., eps=1e-4)
Vfull = SingleLayerPotentialOperator(fesL2, intorder=10, leafsize=fesL2.ndof, eta=3., eps=1e-4)
print ("ndof =", V.mat.nze, "shape =", V.mat.shape)

x = V.mat.CreateRowVector()
x.SetRandom(4711)
err = Norm ( (V.mat-Vfull.mat) * x )
print ("err = ", err)

def test_answer():
    assert V.mat.nze < 5*10**6
    assert err < 1e-6
