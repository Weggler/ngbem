#  testfile for checking integration 
#  output of testrun 29.11.23
# err-integration (8-10):  1.6205859881425085e-07
# err-integration (10-12):  1.5225443842172683e-08
# err-integration (12-14):  1.3555959613067985e-09
# err-integration (14-16):  1.3309124393862712e-10
# err-integration (16-18):  1.507548088767588e-11

# testfile for checking integration accuracy accuracy, exspectation for sp

from netgen.occ import *
from ngsolve import *
from libbem import *


sp = Sphere( (0,0,0), 1)
mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(1)

fesL2 = SurfaceL2(mesh, order=0, dual_mapping=False)
fesH1 = H1(mesh, order=1)
u,v = fesL2.TnT()
uH1, vH1 = fesH1.TnT()

b8 = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=8, method="dense")

b10 = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=10, method="dense")

b12 = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=12, method="dense")

b14 = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=14, method="dense")

b16 = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=16, method="dense")

b18 = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=18, method="dense")


x = b8.mat.CreateRowVector()
x.SetRandom(1)

err1 = Norm ( (b8.mat -b10.mat)*x )
print ("err-integration ( 8-10): ", err1 )
err2 = Norm ( (b10.mat -b12.mat)*x )
print ("err-integration (10-12): ", err2)
err3 = Norm ( (b12.mat -b14.mat)*x )
print ("err-integration (12-14): ", err3)
err4 = Norm ( (b14.mat -b16.mat)*x )
print ("err-integration (14-16): ", err4)
err5 = Norm ( (b16.mat -b18.mat)*x )
print ("err-integration (16-18): ", err5)



def test_answer():
    assert err1 < 1e-6
    assert err2 < 1e-6
    assert err3 < 1e-7
    assert err4 < 1e-8
    assert err5 < 1e-10


