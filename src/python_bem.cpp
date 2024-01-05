#include <solve.hpp>        // everything from ngsolve
#include <fem.hpp>        
#include <python_ngstd.hpp> // python bindigs
#include <cmath>
using namespace ngsolve;

#include "ngbem.hpp"
#include "test_compression.hpp"
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
    if(!method.compare("dense")) {
       leafsize = INT_MAX;
       eta = 0.;
    }
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<LaplaceSLKernel<3>>>(space, space, LaplaceSLKernel<3>(), param);
    
  }, py::arg("space"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
    py::arg("method")="svd", py::arg("testhmatrix")=false);


  m.def("DoubleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space,
                                           int intorder, int leafsize, double eta, double eps,
                                           string method, bool testhmatrix) -> shared_ptr<IntegralOperator<>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<LaplaceDLKernel<3>>>(trial_space, test_space, LaplaceDLKernel<3>(), param);    
  }, py::arg("trial_space"), py::arg("test_space"), py::arg("intorder")=3, py::arg("leafsize")=40,
        py::arg("eta")=2., py::arg("eps")=1e-6,
        py::arg("method")="svd", py::arg("testhmatrix")=false);


  m.def("HypersingularOperator", [](shared_ptr<FESpace> space, int intorder, int leafsize, double eta, double eps,
                                    string method, bool testhmatrix) -> shared_ptr<IntegralOperator<>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<LaplaceSLKernel<3>>>(space, space,
                                                                    make_shared<T_DifferentialOperator<DiffOpBoundaryRot>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpBoundaryRot>>(), 
                                                                    LaplaceSLKernel<3>(), param);
    
  }, py::arg("space"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
        py::arg("method")="svd", py::arg("testhmatrix")=false);
  
  

  m.def("HelmholtzSingleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,
                                                    int intorder, int leafsize, double eta, double eps,
                                                    string method, bool testhmatrix) -> shared_ptr<IntegralOperator<Complex>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<HelmholtzSLKernel<3>>>(trial_space, test_space, HelmholtzSLKernel<3>(kappa), param);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr, py::arg("kappa"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
    py::arg("method")="svd", py::arg("testhmatrix")=false);



  m.def("HelmholtzDoubleLayerPotentialOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,
                                                    int intorder, int leafsize, double eta, double eps,
                                                    string method, bool testhmatrix) -> shared_ptr<IntegralOperator<Complex>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<HelmholtzDLKernel<3>>>(trial_space, test_space, HelmholtzDLKernel<3>(kappa), param);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr, py::arg("kappa"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
    py::arg("method")="svd", py::arg("testhmatrix")=false);


  m.def("HelmholtzCombinedFieldOperator", [](shared_ptr<FESpace> trial_space, shared_ptr<FESpace> test_space, double kappa,
                                             int intorder, int leafsize, double eta, double eps,
                                             string method, bool testhmatrix) -> shared_ptr<IntegralOperator<Complex>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<CombinedFieldKernel<3>>>(trial_space, test_space, CombinedFieldKernel<3>(kappa), param);
    
  }, py::arg("trial_space"), py::arg("test_space")=nullptr, py::arg("kappa"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
    py::arg("method")="svd", py::arg("testhmatrix")=false);


  m.def("MaxwellVecSingleLayerPotentialOperator", [](shared_ptr<FESpace> space, double kappa, 
                                                  int intorder, int leafsize, double eta, double eps,
                                                  string method, bool testhmatrix) -> shared_ptr<IntegralOperator<Complex>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<HelmholtzSLKernel<3>>>(space, space,
                                                                    make_shared<T_DifferentialOperator<DiffOpRotatedTrace>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpRotatedTrace>>(), 
                                                                    HelmholtzSLKernel<3>(kappa), param);
    
  }, py::arg("space"), py::arg("kappa"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
        py::arg("method")="svd", py::arg("testhmatrix")=false);
  
  

  m.def("MaxwellScaSingleLayerPotentialOperator", [](shared_ptr<FESpace> space, double kappa, 
                                                  int intorder, int leafsize, double eta, double eps,
                                                  string method, bool testhmatrix) -> shared_ptr<IntegralOperator<Complex>>
  {
    BEMParameters param({intorder, leafsize, eta, eps, method, testhmatrix});
    return make_unique<GenericIntegralOperator<HelmholtzSLKernel<3>>>(space, space,
                                                                    make_shared<T_DifferentialOperator<DiffOpCurlBoundaryEdge<>>>(),
                                                                    make_shared<T_DifferentialOperator<DiffOpCurlBoundaryEdge<>>>(), 
                                                                    HelmholtzSLKernel<3>(kappa), param);
    
  }, py::arg("space"), py::arg("kappa"), py::arg("intorder")=3, py::arg("leafsize")=40, py::arg("eta")=2., py::arg("eps")=1e-6,
        py::arg("method")="svd", py::arg("testhmatrix")=false);
  
  

  
  
  m.def("TestCompression", [](int nx, int ny, double eta, double eps, string method) -> py::object
  {
    if (method=="svd")
      return py::cast(TestCompressionSVD (nx, ny, eta, eps));
    
    if (method=="tsvd")
      return py::cast(TestCompressionTSVD (nx, ny, eta, eps));

    if (method=="aca")
      return py::cast(TestCompressionACA (nx, ny, eta, eps));
    
    return py::cast(tuple<int,double>{ 0, 0 });
  },
    "methods: 'svd', 'tsvd', to come: 'aca'\n"
    "returns rank, error, time",
    py::arg("nx")=100, py::arg("ny")=100, py::arg("eta")=2.0,
    py::arg("eps")=1e-10, py::arg("method")="svd");
}

