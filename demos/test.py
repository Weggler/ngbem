import sys
sys.path.append('../build')


from ngsolve import *
from netgen.csg import unit_cube
from netgen.geom2d import unit_square

mesh = Mesh(unit_cube.GenerateMesh(maxh=0.5))
print(mesh.GetBoundaries())

fes = SurfaceL2(mesh, order=0)
fes2 = H1(mesh, order=0, definedon=mesh.Boundaries(".*"))

fes2 = Compress (fes2)

a = BilinearForm(fes)
b = BilinearForm(fes2, fes)

u = fes2.TrialFunction()
v = fes.TestFunction()

w = fes2.TestFunction()

from libbem import *

print("pyA")
SingleLayerPotential(a)    # added special element to a
print("pyB")

a.Assemble()
print("pyC")
#print (a.mat)               # 12 x 12


DoubleLayerPotential(b)

print("pyD")
b += 0.5 * u * v * ds
b.Assemble()
print("pyE")
#print (b.mat)               # 12 x 8


gfu = GridFunction(fes2)
gfu.Set (x**2 + y**2 + z**2, definedon=mesh.Boundaries(".*"))
#print (gfu.vec)           

gft = GridFunction(fes)

gft.vec.data += a.mat.Inverse() * (b.mat * gfu.vec)
print(gft.vec)
print(sqrt (Integrate ( gft*gft, mesh, BND)))

