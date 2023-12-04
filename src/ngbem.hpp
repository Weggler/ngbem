#ifndef NGBEM_hpp
#define NGBEM_hpp

#include "hmat.hpp"

namespace ngbem
{
  using namespace ngcomp;

  
  class SingleLayerPotentialOperator : public SpecialElement
  {
    shared_ptr<FESpace> space;
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;
    int intorder;

    ClusterTree cluster_tree;
    Array<Array<int>> elems4dof; // contains list of elems contributing to bnd-dof 


  public:
    SingleLayerPotentialOperator(shared_ptr<FESpace> aspace, int intorder);

    void CalcElementMatrix(FlatMatrix<double> matrix, // matrix dim = ndof_bnd x ndof_bnd
                           LocalHeap &lh) const override;

    void CalcBlockMatrix(FlatMatrix<double> &matrix, Array<DofId> &setI, Array<DofId> &setJ, // matrix dim = size(setI) x size(setJ)
                           LocalHeap &lh);

    void GetDofNrs(Array<int> &dnums) const override;


  };




  class DoubleLayerPotentialOperator : public SpecialElement
  {
    shared_ptr<FESpace> space;
    shared_ptr<FESpace> space2;
    
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;
    Array<DofId> mapglob2bnd2;
    Array<DofId> mapbnd2glob2;
    int intorder;
  public:
    DoubleLayerPotentialOperator(shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace, int intorder);

    void CalcElementMatrix(FlatMatrix<double> matrix, // matrix dim = ndof_bnd x ndof_bnd2
                           LocalHeap &lh) const override;

    void GetDofNrs(Array<int> &dnums) const override; 
    void GetDofNrs2(Array<int> &dnums) const override;
  };

}



#endif

