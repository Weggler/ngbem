# Welcome to NG-BEM

NG-BEM implements boundary integral operators on top of NGSolve.

Currently it supports single-layer and double-layer operators and hypersingular operators for the

* Laplace equation
* Helmholtz equation
* Maxwell equation

NG-BEM supports high order function spaces on curved surface meshes.
It uses numerical integration following [Sauter-Schwab: Boundary Element Methods, 2009], and matrix compression based on the ACA algorithm [Bebendorf: Hierarchical Matrices: A Means to Efficiently Solve Elliptic Boundary Value Problems, 2008]



Installation:

* Install a recent NGSolve (later than Feb. 10, 2024)

* Install ngbem from github:

      git clone https://github.com/Weggler/ngbem.git
      cd ngbem
      mkdir build
      cd build
      cmake ..
      make -j4 install

Try notebooks from the demos folder.
You can read (only read) them directly on github: https://github.com/Weggler/ngbem/tree/main/demos


```{tableofcontents}
```
