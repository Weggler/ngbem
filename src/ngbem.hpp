#ifndef NGBEM_hpp
#define NGBEM_hpp

#include "hmat.hpp"

namespace ngbem
{
  using namespace ngcomp;

  struct BEMParameters
  {
    const int intorder;
    const int leafsize;
    const double eta;
    const double eps;
    const char *method;
  };
  
  class SingleLayerPotentialOperator : public SpecialElement
  {
    shared_ptr<FESpace> space;
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;
    struct BEMParameters param;

    ClusterTree cluster_tree;
    shared_ptr<HMatrix> hmatrix;
    Array<Array<int>> elems4dof; // contains list of elems contributing to global dof 

  public:
    SingleLayerPotentialOperator(shared_ptr<FESpace> aspace, struct BEMParameters param);

    // matrix dim = ndof x ndof, volume dofs are not used
    void CalcElementMatrix(FlatMatrix<double> matrix, 
                           LocalHeap &lh) const override;

    void CalcBlockMatrix(FlatMatrix<double> matrix, const Array<DofId> &trialdofs, const Array<DofId> &testdofs, 
			 LocalHeap &lh) const;

    unique_ptr<LowRankMatrix> CalcFarFieldBlock(const Array<DofId> &trialdofs, const Array<DofId> &testdofs, LocalHeap &lh) const;
    
    void CalcHMatrix(HMatrix & hmatrix, LocalHeap &lh, struct BEMParameters &param) const;

    virtual void Apply (FlatVector<double> elx, FlatVector<double> ely, 
			LocalHeap & lh) const override;
    
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

    ClusterTree cluster_tree;
    Array<Array<int>> elems4dof; // contains list of elems contributing to bnd-dof
    ClusterTree cluster_tree2;
    Array<Array<int>> elems4dof2; // contains list of elems contributing to bnd-dof 

  public:
    DoubleLayerPotentialOperator(shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace, int intorder);

    void CalcElementMatrix(FlatMatrix<double> matrix, // matrix dim = ndof_bnd x ndof_bnd2
                           LocalHeap &lh) const override;
    
    void CalcBlockMatrix(FlatMatrix<double> matrix, const Array<DofId> &testdofs, const Array<DofId> &trialdofs, // matrix dim = size(testdofs) x size(trialdofs)
			 LocalHeap &lh) const;

    void GetDofNrs(Array<int> &dnums) const override; 
    void GetDofNrs2(Array<int> &dnums) const override;
  };

}



#endif

