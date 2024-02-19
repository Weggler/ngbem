#!/usr/bin/env python

import sys
sys.path.append("../build/")
from netgen.occ import *
import netgen.meshing as meshing
from ngsolve import *
from libbem import *
from ngsolve import Projector, Preconditioner
from ngsolve.krylovspace import GMRes
from ngsolve.fem import CompilePythonModule
from pathlib import Path

# Get reference mie series current from cpp
txt = Path('mie.cpp').read_text() 
mie = CompilePythonModule(txt, init_function_name='Mie', add_header=False)
miecurrent = mie.MieCurrent()

# Scattering on a sphere
sp = Sphere((0, 0, 0), 0.25)
kappa = 5.
E_inc = CF((1, 0, 0)) * exp(1j * kappa * z)

print("order ndof error")
for order in range(0, 5):
    for n in range(1, 5):
        h = 0.1 / n
        mesh = Mesh(OCCGeometry(sp).GenerateMesh(maxh=h, perfstepsend=meshing.MeshingStep.MESHSURFACE)).Curve(4)
        fesHDiv = HDivSurface(mesh, order=order, complex=True)
        uHDiv, vHDiv = fesHDiv.TnT() 

        rhs = LinearForm(-E_inc * vHDiv.Trace() * ds(bonus_intorder=10)).Assemble()

        j = GridFunction(fesHDiv)
        pre = BilinearForm(uHDiv.Trace() * vHDiv.Trace() * ds).Assemble().mat.Inverse(freedofs=fesHDiv.FreeDofs()) 
        with TaskManager(): 
            V = MaxwellSingleLayerPotentialOperator(fesHDiv, kappa, intorder=8 + 2 * order, leafsize=320, eta=3., eps=1e-10)
            GMRes(A=V.mat, pre=pre, b=rhs.vec, x=j.vec, tol=1e-11, maxsteps=2000, printrates=False)
    
        j.vec[:] *= kappa
    
        error = sqrt(Integrate(Norm(j - miecurrent)**2, mesh, BND))
        print(order, fesHDiv.ndof, error)
        
