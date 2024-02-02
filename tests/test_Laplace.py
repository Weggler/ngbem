#!/usr/bin/env python
# coding: utf-8

import sys
sys.path.append("../build/")
from ngsolve import *
from netgen.occ import *
from ngsolve.krylovspace import CG, GMRes
from ngsolve.webgui import Draw
from libbem import *

# Reference errors (dirichlet, neumann, screen)
errref = [0.0004, 0.007, 8e-07]

# Bottom sphere: Dirichlet boundary, top sphere: Neumann boundary
topsphere = Sphere((0,0,0), 1) * Box((-1,-1,0),(1,1,1))
botsphere = Sphere((0,0,0), 1) - Box((-1,-1,0),(1,1,1))
topsphere.faces.name = "neumann"
botsphere.faces.name = "dirichlet"
shape = Fuse([topsphere,botsphere])
screen = WorkPlane(Axes( (0,0,0), Z, X)).RectangleC(0.5,0.5).Face()
screen.faces.name="screen"
mesh_screen = Mesh(OCCGeometry(screen).GenerateMesh(maxh=1)).Curve(1)

# h-refinement for fixed order
order = 3

mesh = Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.343)).Curve(order)

fesH1 = H1(mesh, order=order, dirichlet="dirichlet", definedon=mesh.Boundaries(".*"))
u1,v1 = fesH1.TnT()
fesL2 = SurfaceL2(mesh, order=order-1, dirichlet="neumann")
u,v = fesL2.TnT()

uexa = CF(1. / sqrt((x-1)**2 + (y-1)**2 + (z-1)**2))
n = specialcf.normal(3)
gradn_uexa = CF((uexa.Diff(x), uexa.Diff(y), uexa.Diff(z))) * n

# Dirichlet data on Dirichlet boundary
ud = GridFunction(fesH1)
ud.Interpolate(uexa, definedon=mesh.Boundaries("dirichlet"))

# Neumann data on Neumann boundary
un = GridFunction(fesL2)
un.Interpolate(gradn_uexa, definedon=mesh.Boundaries("neumann"))

eps = 1e-9
intorder = 3 * order + 6
with TaskManager():
    V = SingleLayerPotentialOperator(fesL2, intorder=intorder, eps=eps)
    K = DoubleLayerPotentialOperator(fesH1, fesL2, intorder=intorder, eps=eps)
    W = HypersingularOperator(fesH1, intorder=intorder, eps=eps)
    M = BilinearForm(u1.Trace() * v.Trace() * ds(bonus_intorder=3)).Assemble()
    fd = ((0.5 * M.mat + K.mat) * ud.vec - V.mat * un.vec).Evaluate()
    fn = ((0.5 * M.mat.T - K.mat.T) * un.vec - W.mat * ud.vec).Evaluate()
    pred = BilinearForm(u.Trace() * v.Trace() * ds(bonus_intorder=3), check_unused=False).Assemble()
    pren = BilinearForm(u1.Trace() * v1.Trace() * ds(bonus_intorder=3), check_unused=False).Assemble()

    lhs = BlockMatrix([[V.mat, - K.mat], [K.mat.T, W.mat]])
    rhs = BlockVector([fd, fn])
    pre = BlockMatrix ([[pred.mat.Inverse(freedofs=fesL2.FreeDofs()), None], [None, pren.mat.Inverse(freedofs=fesH1.FreeDofs())]])

    sol = GMRes(A=lhs, b=rhs, pre=pre, maxsteps=400, printrates=False, tol=1e-12)

# Neumann data on Dirichlet boundary
gfu1 = GridFunction(fesL2)
gfu1.vec[:] = sol[0]

# Dirichlet data on Neumann boundary
gfu0 = GridFunction(fesH1)
gfu0.vec[:] = sol[1]

# Compute the L2-error in the traces on the complete boundary
errd = sqrt(Integrate((uexa - gfu0 - ud)**2, mesh.Boundaries(".*"), BND))
errn = sqrt(Integrate((gradn_uexa - gfu1 - un)**2, mesh.Boundaries(".*"), BND))

# Post-processing on screen
fes_screen = H1(mesh_screen, order=7)
gf_screen = GridFunction(fes_screen)
with TaskManager():
    gf_screen.Set(V.GetPotential(gfu1) - K.GetPotential(gfu0) + V.GetPotential(un) - K.GetPotential(ud), definedon=mesh_screen.Boundaries("screen"), dual=False)

    errs = sqrt(Integrate((gf_screen - uexa)**2, mesh_screen.Boundaries(".*"), BND))    

def test_answer():
    assert errd < errref[0]
    assert errn < errref[1]
    assert errs < errref[2]
