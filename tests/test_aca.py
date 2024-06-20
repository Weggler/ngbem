from netgen.occ import *
from ngsolve import *
from ngbem import *


sp = Sphere( (0,0,0), 1)
mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(3)

fesL2 = SurfaceL2(mesh, order=3, dual_mapping=True)
u,v = fesL2.TnT()

eps = 1e-6;

with TaskManager(pajetrace=1000*1000*1000):
    Vdense = SingleLayerPotentialOperator(fesL2, intorder=10, method="dense")
    V = SingleLayerPotentialOperator(fesL2, intorder=10, leafsize=40, eta=3., eps=eps, method="aca", )

print (V.mat.nze)
x = V.mat.CreateRowVector()
x.SetRandom(4711)
err = Norm( (Vdense.mat-V.mat)*x )
print ( "hmatrix error =", err)



def test_answer():
    assert V.mat.nze < 15*10**6
    assert err < eps
    
    
