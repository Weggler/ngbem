#include <solve.hpp>        // everything from ngsolve
#include <python_ngstd.hpp> // python bindigs
#include <cmath>
using namespace ngsolve;

#include "ngbem.hpp"
using namespace ngbem;

PYBIND11_MODULE(libbem, m)
{
  cout << "Loading ngbem library" << endl;

  m.def("SingleLayerPotential", [](BilinearForm &bfa, int intorder) {
    // cout << "create single-layer potential" << endl;
    auto bemop = make_unique<SingleLayerPotential>(bfa.GetTrialSpace(), intorder);
    bfa.AddSpecialElement(std::move(bemop));
  }, py::arg("bf"), py::arg("intorder")=10);



  /*
  // was something for testing ???
  m.def("SingleLayerPotential", [](shared_ptr<ProxyFunction> trialfunc,
                                   shared_ptr<ProxyFunction> testfunc) {
    cout << "create single-layer from trial- and testfunction" << endl;

    auto trialspace = trialfunc->GetFESpace();
    auto mesh = trialspace->GetMeshAccess();

    auto evaluator = trialfunc->Evaluator();
    LocalHeap lh(10*1000*1000);
    for (int i = 0; i < mesh->GetNSE(); i++)
      {
        HeapReset hr(lh);
        Array<DofId> dnums;
        trialspace->GetDofNrs(ElementId(BND, i), dnums);
        auto & fel_trial = trialspace->GetFE(ElementId(BND, i), lh);
        IntRange r1 = evaluator->UsedDofs(fel_trial);
        cout << "dofs = " << dnums << endl;
        cout << "subdofs = " << dnums[r1] << endl;
        
        // for (auto d : dnums) // modern C++11 loop
        // bnddofs.SetBit(d);
      }

    // auto bemop = make_unique<BEMOperator>(bfa.GetTrialSpace());
    // bfa.AddSpecialElement(move(bemop));
  }, py::arg("trialfunc"), py::arg("testfunc"));
  */


  m.def("DoubleLayerPotential", [](BilinearForm &bfb, int intorder) {
    // cout << "create double-layer potential" << endl;

    auto bemop2 = make_unique<DoubleLayerPotential>(bfb.GetTrialSpace(), bfb.GetTestSpace(), intorder);
    bfb.AddSpecialElement(std::move(bemop2));
  }, py::arg("bf"), py::arg("intorder")=10);
  ;
}

