#include <solve.hpp>        // everything from ngsolve
#include <python_ngstd.hpp> // python bindigs
#include <cmath>
using namespace ngsolve;

#include "ngbem.hpp"
using namespace ngbem;

PYBIND11_MODULE(libbem, m)
{
  cout << "Loading ngbem library" << endl;

  m.def("SingleLayerPotentialOperator", [](BilinearForm &bfa, int intorder, int leafsize, double eta, double eps, char *method) {
    // cout << "create single-layer potential" << endl;
    struct BEMParameters param({intorder, leafsize, eta, eps, method});
    auto bemop = make_unique<SingleLayerPotentialOperator>(bfa.GetTrialSpace(), param);
    bfa.AddSpecialElement(std::move(bemop));
  }, py::arg("bf"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6, py::arg("method")="svd");


  m.def("DoubleLayerPotentialOperator", [](BilinearForm &bfb, int intorder) {
    // cout << "create double-layer potential" << endl;

    auto bemop2 = make_unique<DoubleLayerPotentialOperator>(bfb.GetTrialSpace(), bfb.GetTestSpace(), intorder);
    bfb.AddSpecialElement(std::move(bemop2));
  }, py::arg("bf"), py::arg("intorder")=10);
  ;
}

