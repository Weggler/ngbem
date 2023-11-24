#include <solve.hpp>        // everything from ngsolve
#include <python_ngstd.hpp> // python bindigs
#include <cmath>
using namespace ngsolve;

#include "ngbem.hpp"
using namespace ngbem;

PYBIND11_MODULE(libbem, m)
{
  cout << "Loading ngbem library" << endl;

  m.def("SingleLayerPotentialOperator", [](BilinearForm &bfa, int intorder) {
    // cout << "create single-layer potential" << endl;
    auto bemop = make_unique<SingleLayerPotentialOperator>(bfa.GetTrialSpace(), intorder);
    bfa.AddSpecialElement(std::move(bemop));
  }, py::arg("bf"), py::arg("intorder")=10);


  m.def("DoubleLayerPotentialOperator", [](BilinearForm &bfb, int intorder) {
    // cout << "create double-layer potential" << endl;

    auto bemop2 = make_unique<DoubleLayerPotentialOperator>(bfb.GetTrialSpace(), bfb.GetTestSpace(), intorder);
    bfb.AddSpecialElement(std::move(bemop2));
  }, py::arg("bf"), py::arg("intorder")=10);
  ;
}

