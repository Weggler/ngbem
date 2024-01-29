# testfile for checking integration accuracy accuracy, exspectation for dp

import sys
sys.path.append("../build/")
from netgen.occ import *
from ngsolve import *
from libbem import *

errref = [0.00014356401693177312,
          1.1556523998207007e-05,
          9.959433217687845e-07,
          9.474259269075394e-08,
          9.6163481617507e-09,
          9.182994523208552e-10,
          8.329321209219403e-11,
          7.326096658552548e-12,
          8.665555899250029e-13,
          1.6241027931801214e-13]

def power_iteration(A, x, y):
    z = x
    for i in range(30):
        x = (1. / Norm(z)) * z
        y = A * x
        z = A.T * y
    return y;

def test_answer():
    sp = Sphere( (0,0,0), 1)
    mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=0.5)).Curve(1)
    fesL2 = SurfaceL2(mesh, order=0, dual_mapping=False)
    fesH1 = H1(mesh, order=1)
    aref = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=26, method="dense")
    x = aref.mat.CreateRowVector()
    y = aref.mat.CreateColVector()

    for n in range(10):
        intorder = 2 * (n + 1)
        a = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=intorder, method="dense")
        x.SetRandom(1)
        z = power_iteration(aref.mat - a.mat, x, y)
        assert Norm(z ) < 5. * errref[n]

