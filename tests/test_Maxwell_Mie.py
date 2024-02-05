#!/usr/bin/env python

import sys
sys.path.append("../build/")
from netgen.occ import *
import netgen.meshing as meshing
import netgen.gui
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
order = 4
sp = Sphere((0, 0, 0), 0.25)
mesh = Mesh(OCCGeometry(sp).GenerateMesh(maxh=1, perfstepsend=meshing.MeshingStep.MESHSURFACE)).Curve(order)

fesHDiv = HDivSurface(mesh, order=order, complex=True)
uHDiv,vHDiv = fesHDiv.TnT() # H(div_Gamma) trial space for Neumann data ( nx curlE ) and test space for BIE
print ("ndof HDiv =", fesHDiv.ndof)

kappa = 5.
E_inc = CF((1, 0, 0)) * exp(1j * kappa * z)
rhs = LinearForm(-E_inc * vHDiv.Trace() * ds(bonus_intorder=10)).Assemble()  # <H(curl_Gamma), H(div_Gamma)> 

j = GridFunction(fesHDiv)
pre = BilinearForm(uHDiv.Trace() * vHDiv.Trace() * ds).Assemble().mat.Inverse(freedofs=fesHDiv.FreeDofs()) 
with TaskManager(): 
    V = MaxwellSingleLayerPotentialOperator(fesHDiv, kappa, intorder=16, leafsize=320, eta=0., eps=1e-8)
    GMRes(A=V.mat, pre=pre, b=rhs.vec, x=j.vec, tol=1e-8, maxsteps=1000, printrates=True)
    
j.vec[:] *= kappa
    
error = sqrt(Integrate(Norm(j - miecurrent)**2, mesh, BND))
print("L2-error in j: ", error)

Draw(miecurrent, mesh, "mie", draw_vol=False, order=order)
Draw(j, mesh, "j", draw_vol=False, order=order)
Draw(Norm(j - miecurrent), mesh, "j - mie", draw_vol=False, order=order)

jex = GridFunction(fesHDiv)
jex.Set(miecurrent)
print("jex-norm", Norm(jex.vec))
jex.vec[:] *= 1./ kappa
res = (V.mat * jex.vec).Evaluate() - rhs.vec
print("residual: ", Norm(res) / Norm(rhs.vec))
