#!/usr/bin/env python
# coding: utf-8

from ngsolve import *
from netgen.occ import *
from ngsolve.krylovspace import GMRes
from ngbem import *

# reference errors vol, dirichlet, neumann for costabel and nedelec coupling
errref = [0.0006,
          0.0007,
          0.002,
          0.0006,
          0.0007,
          0.002]

order = 2
intorder = 16

radius = 1
shape = Sphere((0,0,0), radius)
mesh = Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.5)).Curve(order)

fesH1 = H1(mesh, order=order)
u,v = fesH1.TnT()
a = BilinearForm((grad(u)*grad(v))*dx).Assemble()
source = 3. / (4. * pi * radius ** 3)
f = LinearForm(source*v*dx()).Assemble()

fesL2 = SurfaceL2(mesh, order=order-1)
f2 = LinearForm(fesL2).Assemble()  # 0-vector

with TaskManager():
    V = SingleLayerPotentialOperator(fesL2, intorder=intorder, eta=0.)
    K = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=intorder, eta=0.)
    D = HypersingularOperator(fesH1, intorder=intorder, eta=0.)
    M = BilinearForm(fesH1.TrialFunction()*fesL2.TestFunction().Trace()*ds(bonus_intorder=10)).Assemble()

costabel = BlockMatrix ( [ [ a.mat+D.mat, (-0.5*M.mat+K.mat).T ], [ (-0.5*M.mat+K.mat), -V.mat ] ] )
nedelec = BlockMatrix ( [ [ a.mat, -M.mat.T ], [ (-0.5*M.mat+K.mat), -V.mat ] ] )
rhs = BlockVector( [ f.vec, f2.vec ] )

bndmass = BilinearForm( fesL2.TrialFunction().Trace()*fesL2.TestFunction().Trace()*ds(bonus_intorder=10)).Assemble()
h1mass = BilinearForm( fesH1.TrialFunction().Trace()*fesH1.TestFunction().Trace()*dx).Assemble()
pre = BlockMatrix ( [ [h1mass.mat.Inverse(), None], [None, bndmass.mat.Inverse(freedofs=fesL2.FreeDofs())] ])

with TaskManager():
    sol_costabel = GMRes(A=costabel, b=rhs, pre=pre, tol=1e-10, maxsteps=400, printrates=False)
    sol_nedelec = GMRes(A=nedelec, b=rhs, pre=pre, tol=1e-10, maxsteps=400, printrates=False)

uexa = 1. / (4. * pi) * (3. / (2. * radius) - (x**2 + y**2 + z**2) / (2. * radius ** 3))
n = specialcf.normal(3)
gradn_uexa = CF((uexa.Diff(x), uexa.Diff(y), uexa.Diff(z))) * n
u0 = GridFunction(fesH1)
u0.Interpolate(uexa)
u1 = GridFunction(fesL2)
u1.Interpolate(gradn_uexa, definedon=mesh.Boundaries(".*"))

def test_answer():
    gfu0 = GridFunction(fesH1)
    gfu1 = GridFunction(fesL2)
    gfu0.vec[:] = sol_costabel[0]
    gfu1.vec[:] = sol_costabel[1]
    assert sqrt(Integrate((gfu0 - u0) ** 2, mesh)) < errref[0]
    assert sqrt(Integrate((gfu0 - u0) ** 2, mesh.Boundaries(".*"), BND)) < errref[1]
    assert sqrt(Integrate((gfu1 - u1) ** 2, mesh.Boundaries(".*"), BND)) < errref[2]

    gfu0.vec[:] = sol_nedelec[0]
    gfu1.vec[:] = sol_nedelec[1]
    assert sqrt(Integrate((gfu0 - u0) ** 2, mesh)) < errref[3]
    assert sqrt(Integrate((gfu0 - u0) ** 2, mesh.Boundaries(".*"), BND)) < errref[4]
    assert sqrt(Integrate((gfu1 - u1) ** 2, mesh.Boundaries(".*"), BND)) < errref[5]

