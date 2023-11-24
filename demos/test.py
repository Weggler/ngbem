# Differentialoperator: Laplace Operator  
# Boundary conditions: Dirichlet boundary conditions
# Domain: bounded, 3D 
# Solving for Neumann data using BEM 

import sys
sys.path.append('../build')

from ngsolve import *
from netgen.csg import unit_cube
from netgen.geom2d import unit_square

from libbem import *

# define the mesh
mesh = Mesh(unit_cube.GenerateMesh(maxh=0.5))
#print(mesh.GetBoundaries())

# define finite element spaces
fes = SurfaceL2(mesh, order=0)
fes2 = H1(mesh, order=0, definedon=mesh.Boundaries(".*"))
fes2 = Compress (fes2)
u = fes2.TrialFunction()
v = fes.TestFunction()
w = fes2.TestFunction()

# assemble single layer potential operator matrix
a = BilinearForm(fes)
SingleLayerPotentialOperator(a, intorder = 12)  # added special element to a
a.Assemble()

# assemble double layer potential operator matrix
b = BilinearForm(fes2, fes)
DoubleLayerPotentialOperator(b, intorder = 12)
b += - 0.5 * u * v.Trace() * ds(bonus_intorder=3)
b.Assemble()

# define dirichlet data
gfu = GridFunction(fes2)
gfu.Set (x**2 + y**2 + z**2, definedon=mesh.Boundaries(".*"))

# compute neumann data
gft = GridFunction(fes)
gft.vec.data = a.mat.Inverse() * (b.mat * gfu.vec)
#print(gft.vec)
#print(sqrt (Integrate ( gft*gft, mesh, BND)))

