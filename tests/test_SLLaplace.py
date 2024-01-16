from netgen.occ import *
from ngsolve import *
from libbem import *


sp = Sphere( (0,0,0), 1)
mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(3)

fesL2 = SurfaceL2(mesh, order=3, dual_mapping=True)
u,v = fesL2.TnT()

with TaskManager(pajetrace=1000*1000*1000):
    Vdense = SingleLayerPotentialOperator(fesL2, intorder=10, method="dense")
    V = SingleLayerPotentialOperator(fesL2, intorder=10, method="svd")

print (V.mat.nze)
x = V.mat.CreateRowVector()
x.SetRandom(4711)
err = Norm( (Vdense.mat-V.mat)*x )
print ( "hmatrix error =", err)



def test_answer():
    assert V.mat.nze < 15*10**6
    assert err < 1e-8
    
    
