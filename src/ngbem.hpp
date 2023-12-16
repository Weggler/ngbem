#ifndef NGBEM_hpp
#define NGBEM_hpp

#include "hmat.hpp"

namespace ngbem
{
  using namespace ngcomp;

  /** BEMParameters. 
      All parameters to define approximation of the layer potential operators. */
  struct BEMParameters
  {
    /* intorder defines integration rule for singular pairings. */  
    const int intorder;
    /* leafsize is the minimal size of a #Custer */
    const int leafsize;
    /* eta defines admissibily condition for cluster pairs*/
    const double eta;
    /* eps defines low rank approximation */
    const double eps;
    /* method defines the low rank approximation method */
    const string method;
    // testing hmatrix accuracy
    bool testhmatrix;
  };
  
  /** SingleLayerPotentialOperator. 
      */
  class SingleLayerPotentialOperator : public SpecialElement
  {
    shared_ptr<FESpace> space; 
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;
    Array<Array<int>> elems4dof; // contains list of elems contributing to dof 

    BEMParameters param;
    shared_ptr<ClusterTree> cluster_tree;
    shared_ptr<HMatrix> hmatrix;

  public:
    /** Constructor */
    SingleLayerPotentialOperator(shared_ptr<FESpace> aspace, struct BEMParameters param);

    /** Routine computes the SL potential matrix according to the given FE space. 
        The matrix is dense with dim = ndof(L2) x ndof(L2), where volume dofs are not used */
    void CalcElementMatrix(FlatMatrix<double> matrix, LocalHeap &lh) const override;

    /** Compute the sub-block of SL potential matrix which belongs to the given dofs, 
        where #matrix is dense. We use the routine to compute a #BEMBlock of type "nearfield". */
    void CalcBlockMatrix(FlatMatrix<double> matrix, const Array<DofId> &trialdofs, const Array<DofId> &testdofs, 
			 LocalHeap &lh) const;

    /** Compute the sub-block of SL potential matrix which belongs to the given dofs,  
        where #matrix is a #LowRankMatrix. We use the routine to compute a #BEMBlock of type "farfield". */
    unique_ptr<LowRankMatrix> CalcFarFieldBlock(const Array<DofId> &trialdofs, const Array<DofId> &testdofs, LocalHeap &lh) const;
    
    /** Given a #HMatrix structure, compute all block. Memory for farfield blocks gets clear here. */
    void CalcHMatrix(HMatrix & hmatrix, LocalHeap &lh, struct BEMParameters &param) const;

    /** Compute the matrix-vector multiplication $y = A x$. Apply is used by the iterative solver. */
    virtual void Apply (FlatVector<double> elx, FlatVector<double> ely, LocalHeap & lh) const override;
    
    /** Get list of boundary degrees of freedom of given FE space.  */
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

    struct BEMParameters param;

    shared_ptr<ClusterTree> cluster_tree;
    Array<Array<int>> elems4dof; // contains list of elems contributing to bnd-dof
    shared_ptr<ClusterTree> cluster_tree2;
    Array<Array<int>> elems4dof2; // contains list of elems contributing to bnd-dof 
    shared_ptr<HMatrix> hmatrix;

  public:
    DoubleLayerPotentialOperator(shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace, struct BEMParameters param);

    /** Routine computes the DL potential matrix according to the given FE space. 
        The matrix is dense with dim = ndof(L2) x ndof(H1), where volume dofs are not used */
    void CalcElementMatrix(FlatMatrix<double> matrix, LocalHeap &lh) const override;
    
    /** Compute the sub-block of DL potential matrix which belongs to the given dofs, 
        where #matrix is dense with dim = size(testdofs) x size(trialdofs). 
        We use the routine to compute a #BEMBlock of type "nearfield". */
    void CalcBlockMatrix(FlatMatrix<double> matrix, const Array<DofId> &trialdofs, const Array<DofId> &testdofs, 
			 LocalHeap &lh) const;

    /** Compute the sub-block of SL potential matrix which belongs to the given dofs,  
        where #matrix is a #LowRankMatrix. We use the routine to compute a #BEMBlock of type "farfield". */
    unique_ptr<LowRankMatrix> CalcFarFieldBlock(const Array<DofId> &trialdofs, const Array<DofId> &testdofs, LocalHeap &lh) const;
    
    /** Given a #HMatrix structure, compute all block. Memory for farfield blocks gets clear here. */
    void CalcHMatrix(HMatrix & hmatrix, LocalHeap &lh, struct BEMParameters &param) const;

    /** Compute the matrix-vector multiplication $y = A x$. Apply is used by the iterative solver. */
    virtual void Apply (FlatVector<double> elx, FlatVector<double> ely, LocalHeap & lh) const override;

    /** Get list of boundary degrees of freedom of trial space.  */
    void GetDofNrs(Array<int> &dnums) const override; 

    /** Get list of boundary degrees of freedom of test space.  */
    void GetDofNrs2(Array<int> &dnums) const override;
  };

}



#endif

