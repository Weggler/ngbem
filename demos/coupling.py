import sys
sys.path.append('../build')


import numpy as np
import scipy.sparse as sp
from ngsolve import *
from netgen.csg import unit_cube
from netgen.geom2d import unit_square
from libbem import *


mesh = Mesh(unit_cube.GenerateMesh(maxh=0.1))

surfL2 = SurfaceL2(mesh, order=0)
fesH1 = H1(mesh, order=0)

V = BilinearForm(surfL2)
KI = BilinearForm(fesH1, surfL2)
A = BilinearForm(fesH1)
B = BilinearForm(surfL2, fesH1)

rhs1 = LinearForm(fesH1)
rhs2 = LinearForm(surfL2)

sol1 = GridFunction(fesH1)
sol2 = GridFunction(surfL2)

us = surfL2.TrialFunction()
vs = surfL2.TestFunction()

uH1 = fesH1.TrialFunction()
vH1 = fesH1.TestFunction()

SingleLayerPotential(V)
DoubleLayerPotential(KI)

A += grad(uH1) * grad(vH1) * dx 
B += - us * vH1 * ds
KI += - 0.5 * uH1 * vs * ds

rhs1 += vH1 * dx
# rhs = 0

A.Assemble()
B.Assemble()
KI.Assemble()
V.Assemble()

rhs1.Assemble()
rhs2.Assemble()

rows,cols,vals = A.mat.COO()
np_A = sp.csr_matrix((vals,(rows,cols)), shape=A.mat.shape).todense()

rows,cols,vals = B.mat.COO()
np_B = sp.csr_matrix((vals,(rows,cols)), shape=B.mat.shape).todense()

rows,cols,vals = KI.mat.COO()
np_KI = sp.csr_matrix((vals,(rows,cols)), shape=KI.mat.shape).todense()

rows,cols,vals = V.mat.COO()
np_V = sp.csr_matrix((vals,(rows,cols)), shape=V.mat.shape).todense()

np_rhs1=rhs1.vec.FV().NumPy()
np_rhs2=rhs2.vec.FV().NumPy()

print(A.mat.shape)
print(np_A.shape)
print(B.mat.shape)
print(np_B.shape)
print(KI.mat.shape)
print(np_KI.shape)
print(V.mat.shape)
print(np_V.shape)

M = np.block([[np_A,np_B], [-np_KI, np_V]])
rhs = np.block([np_rhs1, np_rhs2])
sol = np.linalg.solve(M, rhs)

sol1.vec.FV().NumPy()[:] = sol[0:np_rhs1.shape[0]]
sol2.vec.FV().NumPy()[:] = sol[np_rhs1.shape[0]:]

print(sol1.vec)
print(sol2.vec)

print(sqrt (Integrate ( sol1 * sol1, mesh)))
print(sqrt (Integrate ( sol2 * sol2, mesh, BND)))



# we want to solve a system [ A  B] x [u] = [rhs]
#                           [-KI V] x [t] = [0  ]
