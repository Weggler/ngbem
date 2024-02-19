#!/usr/bin/env python

import sys
sys.path.append("../build/")
from netgen.occ import *
import netgen.meshing as meshing
from ngsolve import *
from libbem import *
from ngsolve.krylovspace import GMRes
from ngsolve.fem import CompilePythonModule
from pathlib import Path

# reference L2-error in current
errref = 0.0015

# Get reference mie series current from cpp
txt = Path('mie.cpp').read_text() 
mie = CompilePythonModule(txt, init_function_name='Mie', add_header=False)
miecurrent = mie.MieCurrent()

# Scattering on a sphere
order = 4
sp = Sphere((0, 0, 0), 0.25)
mesh = Mesh(OCCGeometry(sp).GenerateMesh(maxh=1, perfstepsend=meshing.MeshingStep.MESHSURFACE)).Curve(order)

fesHDiv = HDivSurface(mesh, order=order, complex=True)
uHDiv,vHDiv = fesHDiv.TnT() # H(div_Gamma) trial space for Neumann data ( nx curlE ) and test space for BIE

kappa = 5.
E_inc = CF((1, 0, 0)) * exp(1j * kappa * z)
rhs = LinearForm(-E_inc * vHDiv.Trace() * ds(bonus_intorder=10)).Assemble()  # <H(curl_Gamma), H(div_Gamma)> 

j = GridFunction(fesHDiv)
pre = BilinearForm(uHDiv.Trace() * vHDiv.Trace() * ds).Assemble().mat.Inverse(freedofs=fesHDiv.FreeDofs()) 
with TaskManager(): 
    V = MaxwellSingleLayerPotentialOperator(fesHDiv, kappa, intorder=16, leafsize=120, eta=3., eps=1e-8)
    GMRes(A=V.mat, pre=pre, b=rhs.vec, x=j.vec, tol=1e-8, maxsteps=1000, printrates=False)
    
j.vec[:] *= kappa
    
error = sqrt(Integrate(Norm(j - miecurrent)**2, mesh, BND))

def test_answer():
    assert error < errref
