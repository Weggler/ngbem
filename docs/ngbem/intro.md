# Welcome to NG-BEM

NG-BEM implements boundary integral operators on top of NGSolve.

Currently it supports single-layer and double-layer operators for the

* Laplace equation
* Helmholtz equation
* Maxwell equation

It uses matrix compression based on the ACA algorithm.



Installation:

* Install a recent NGSolve (later than Jan. 5, 2024)

* Install ngbem from github:

      git clone https://github.com/Weggler/ngbem.git
      cd ngbem
      mkdir build
      cd build
      cmake ..
      make -j4 install

Try notebooks from the demos folder.

```{tableofcontents}
```
