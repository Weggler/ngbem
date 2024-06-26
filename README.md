ngbem: Add-on to NGSolve implementing boundary element methods 

Online documentation: https://weggler.github.io/ngbem/intro.html

ngbem is an add-on library to Netgen/NGSolve, a general purpose, high performance finite element library for the numerical solution of partial differential equations. 

The add-on enables the used of the higher order element methods. Both, Netgen/NGSolve and ngbem are written in C++/C with a small Python interface through which it is typically used. ngbem is work an academic software. Goal of the project is to offer fast boundary element methods for all standard problems. 

The current status comprises an Python interface to solve boundary value problems for homogeneous Laplace, Helmholtz and Maxwell equations.


Quick installation: (skip first line in case ngsolve is already installed)

    pip3 install --pre ngsolve
    pip3 install cmake scikit-build-core pybind11_stubgen numpy webgui_jupyter_widgets
    python3 -m pip install --no-build-isolation git+https://github.com/Weggler/ngbem.git
