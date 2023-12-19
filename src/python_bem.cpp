#include <solve.hpp>        // everything from ngsolve
#include <python_ngstd.hpp> // python bindigs
#include <cmath>
using namespace ngsolve;

#include "ngbem.hpp"
using namespace ngbem;

PYBIND11_MODULE(libbem, m)
{
  cout << "Loading ngbem library" << endl;


  py::class_<IntegralOperator<double>, shared_ptr<IntegralOperator<double>>> (m, "IntegralOperator")
    .def_property_readonly("mat", &IntegralOperator<double>::GetMatrix)
    ;
  py::class_<IntegralOperator<Complex>, shared_ptr<IntegralOperator<Complex>>> (m, "IntegralOperatorC")
    .def_property_readonly("mat", &IntegralOperator<Complex>::GetMatrix)
    ;


  
  m.def("SingleLayerPotentialOperator", [](shared_ptr<FESpace> space, int intorder, int leafsize, double eta, double eps,
                                           string method, bool testhmatrix) -> shared_ptr<IntegralOperator<>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    // return make_unique<SingleLayerPotentialOperator>(space, param);

    return make_unique<GenericIntegralOperator<LaplaceSLKernel<3>>>(space, space, LaplaceSLKernel<3>(), param);
    
  }, py::arg("space"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
    py::arg("method")="svd", py::arg("testhmatrix")=false);


  m.def("DoubleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space,
                                           int intorder, int leafsize, double eta, double eps,
                                           string method, bool testhmatrix) -> shared_ptr<IntegralOperator<>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    // return make_shared<DoubleLayerPotentialOperator>(trial_space, test_space, param);
    return make_unique<GenericIntegralOperator<LaplaceDLKernel<3>>>(trial_space, test_space, LaplaceDLKernel<3>(), param);    
  }, py::arg("trial_space"), py::arg("test_space"), py::arg("intorder")=3, py::arg("leafsize")=40,
        py::arg("eta")=2., py::arg("eps")=1e-6,
        py::arg("method")="svd", py::arg("testhmatrix")=false);



  m.def("HelmholtzSingleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,
                                                    int intorder, int leafsize, double eta, double eps,
                                                    string method, bool testhmatrix) -> shared_ptr<IntegralOperator<Complex>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<HelmholtzSLKernel<3>>>(trial_space, test_space, HelmholtzSLKernel<3>(kappa), param);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr, py::arg("kappa"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
    py::arg("method")="svd", py::arg("testhmatrix")=false);

  
}

