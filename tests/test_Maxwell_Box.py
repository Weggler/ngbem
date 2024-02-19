#!/usr/bin/env python

import sys
sys.path.append("../build/")
from netgen.occ import *
import netgen.meshing as meshing
from ngsolve import *
from ngsolve.krylovspace import GMRes
from libbem import *

# reference L2-error in current
errref = 0.015

# Scattering on a cube with manufactured solution
order = 3
box = Box((-5., -1., -1.), (5., 1., 1.)) + Box((-4., -0.5, 1.), (-1., 0.5, 2.)) + Box((3., -0.3, 1.), (4., 0.3, 3.))
#box = Box((-5., -1., -1.), (5., 1., 1.))
h = 1.
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=h, perfstepsend=meshing.MeshingStep.MESHSURFACE))

fesHCurl = HCurl(mesh, order=order, complex=True)
uHCurl, vHCurl = fesHCurl.TnT()
fesHDiv = HDivSurface(mesh, order=order, complex=True)
uHDiv, vHDiv = fesHDiv.TnT()

# Manufactured solution (kernel with source in s and in direction of e)
kappa = 1. / h
s = CF((0., 0., 0.))
e = (1.192, -0.189, 2.745)
e = CF(e) / sqrt(e[0] ** 2 + e[1] ** 2 + e[2] ** 2)
xms = CF((x, y , z)) - s
r = Norm(xms)
f = (-kappa ** 2 * r ** 2 - 1j * kappa * r + 1.) / r ** 3
g = (e * xms) * (-kappa ** 2 * r ** 2 - 3. * 1j * kappa * r + 3.) / r ** 5
E = exp(1j * kappa * r) * (f * e - g * xms)
curlE = CF((E[2].Diff(y) - E[1].Diff(z), E[0].Diff(z) - E[2].Diff(x), E[1].Diff(x) - E[0].Diff(y)))

# Dirichlet data
n = specialcf.normal(3)
mex = -Cross(Cross(n, E), n)
m = GridFunction(fesHCurl)
m.Set(mex, definedon=mesh.Boundaries(".*"), dual=True)
error = sqrt(Integrate(Norm(mex - m) ** 2, mesh, BND) / Integrate(Norm(mex) ** 2, mesh, BND))

j = GridFunction(fesHDiv)
pre = BilinearForm(uHDiv.Trace() * vHDiv.Trace() * ds).Assemble().mat.Inverse(freedofs=fesHDiv.FreeDofs()) 
with TaskManager(): 
    V = MaxwellSingleLayerPotentialOperator(fesHDiv, kappa, intorder=16, leafsize=320, eta=0., eps=1e-8)
    K = MaxwellDoubleLayerPotentialOperator(fesHCurl, fesHDiv, kappa, intorder=16, leafsize=320, eta=0., eps=1e-8)
    M = BilinearForm(uHCurl.Trace() * vHDiv.Trace()* ds(bonus_intorder=3)).Assemble()
    rhs = ((-0.5 * M.mat + K.mat) * m.vec).Evaluate() 
    GMRes(A=V.mat, pre=pre, b=rhs, x=j.vec, tol=1e-8, maxsteps=1000, printrates=False)

# Compute error with respect to exact current
jex = Cross(n, 1/kappa * curlE)
error = sqrt(Integrate(Norm(j - jex) ** 2, mesh, BND) / Integrate(Norm(jex) ** 2, mesh, BND))

def test_answer():
    assert error < errref
