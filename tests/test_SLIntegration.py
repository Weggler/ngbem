# testfile for checking integration accuracy accuracy, exspectation for sp

from netgen.occ import *
from ngsolve import *
from libbem import *


sp = Sphere( (0,0,0), 1)
mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(1)

fes = SurfaceL2(mesh, order=3, dual_mapping=True)
u,v = fes.TnT()

a8 = SingleLayerPotentialOperator(fes, intorder=8, method="dense")

a10 = SingleLayerPotentialOperator(fes, intorder=10, method="dense")

a12 = SingleLayerPotentialOperator(fes, intorder=12, method="dense")

a14 = SingleLayerPotentialOperator(fes, intorder=14, method="dense")

a16 = SingleLayerPotentialOperator(fes, intorder=16, method="dense")

a18 = SingleLayerPotentialOperator(fes, intorder=18, method="dense")


x = a8.mat.CreateRowVector()
x.SetRandom(1)

err1 = Norm ( (a8.mat -a10.mat)*x )
print ("err-integration ( 8-10): ", err1 )
err2 = Norm ( (a10.mat -a12.mat)*x )
print ("err-integration (10-12): ", err2)
err3 = Norm ( (a12.mat -a14.mat)*x )
print ("err-integration (12-14): ", err3)
err4 = Norm ( (a14.mat -a16.mat)*x )
print ("err-integration (14-16): ", err4)
err5 = Norm ( (a16.mat -a18.mat)*x )
print ("err-integration (16-18): ", err5)


def test_answer():
    assert err1 < 1e-5
    assert err2 < 1e-6
    assert err3 < 1e-7
    assert err4 < 1e-8
    assert err5 < 1e-9

