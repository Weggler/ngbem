# testfile for checking integration accuracy accuracy, exspectation for sp

import sys
sys.path.append("../build/")
from netgen.occ import *
from ngsolve import *
from libbem import *

errref = [1.8017794920533596,
          0.02077674341292797,
          0.007395887818216881,
          8.420761045817377e-06,
          1.2916106824050212e-07,
          1.495777703154113e-08,
          1.7312930021482235e-09,
          1.9893964711682858e-10,
          2.2615647594836697e-11,
          2.6809401671932852e-12]

# power iteration to estimate the error in the spectral norm
def power_iteration(A, x):
    y = x
    for i in range(30):
        x = (1. / Norm(y)) * y
        y = A * x
    return y;

def test_answer():
    sp = Sphere( (0,0,0), 1)
    mesh = Mesh( OCCGeometry(sp).GenerateMesh(maxh=1.)).Curve(2)
    fes = SurfaceL2(mesh, order=3, dual_mapping=True)
    aref = SingleLayerPotentialOperator(fes, intorder=26, method="dense")
    x = aref.mat.CreateRowVector()

    # assemble slp with increasing integration order and compute the 2-norm error
    for n in range(10):
        intorder = 2 * (n + 1)
        a = SingleLayerPotentialOperator(fes, intorder=intorder, method="dense")
        x.SetRandom(1)
        y = power_iteration(aref.mat - a.mat, x)
        assert Norm(y) < 5. * errref[n]
        
