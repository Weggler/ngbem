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



  class IntegralOperator : public SpecialElement
  {
  protected:
    shared_ptr<FESpace> trial_space; 
    shared_ptr<FESpace> test_space; 
    
    BEMParameters param;

    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;
    Array<DofId> mapglob2bnd2;
    Array<DofId> mapbnd2glob2;

    Table<int> elems4dof; // contains list of elems contributing to bnd-dof
    Table<int> elems4dof2; // contains list of elems contributing to bnd-dof
    
    shared_ptr<ClusterTree> trial_ct; 
    shared_ptr<ClusterTree> test_ct;
    
    shared_ptr<HMatrix> hmatrix;


  public:
    IntegralOperator (shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space, BEMParameters param);
    shared_ptr<BaseMatrix> GetMatrix() const { return hmatrix; }

    /** Given a #HMatrix structure, compute all block. Memory for farfield blocks gets clear here. */
    void CalcHMatrix(HMatrix & hmatrix, LocalHeap &lh, struct BEMParameters &param) const;

    virtual void CalcBlockMatrix(FlatMatrix<double> matrix,
                                 FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
                                 LocalHeap &lh) const = 0;
    
    /** Compute the sub-block of integral operator matrix which belongs to the given dofs,  
        where #matrix is a #LowRankMatrix. We use the routine to compute a #BEMBlock of type "farfield". */
    virtual unique_ptr<LowRankMatrix>
    CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                      LocalHeap &lh) const;


    /** Routine computes the DL potential matrix according to the given FE space. 
        The matrix is dense with dim = ndof(L2) x ndof(H1), where volume dofs are not used */
    void CalcElementMatrix(FlatMatrix<double> matrix, LocalHeap &lh) const override;


    /** Compute the matrix-vector multiplication $y = A x$. Apply is used by the iterative solver. */
    virtual void Apply (FlatVector<double> elx, FlatVector<double> ely, LocalHeap & lh) const override;
    
    /** Get list of boundary degrees of freedom of trial space.  */
    void GetDofNrs(Array<int> &dnums) const override; 

    /** Get list of boundary degrees of freedom of test space.  */
    void GetDofNrs2(Array<int> &dnums) const override;
  };

  
  /** SingleLayerPotentialOperator. 
      */
  class SingleLayerPotentialOperator : public IntegralOperator
  {
  public:
    /** Constructor */
    SingleLayerPotentialOperator(shared_ptr<FESpace> aspace, BEMParameters param);

    /** Compute the sub-block of SL potential matrix which belongs to the given dofs, 
        where #matrix is dense. We use the routine to compute a #BEMBlock of type "nearfield". */
    void CalcBlockMatrix(FlatMatrix<double> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
			 LocalHeap &lh) const override;
  };


  class DoubleLayerPotentialOperator : public IntegralOperator
  {
  public:
    DoubleLayerPotentialOperator(shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace, struct BEMParameters param);

    
    /** Compute the sub-block of DL potential matrix which belongs to the given dofs, 
        where #matrix is dense with dim = size(testdofs) x size(trialdofs). 
        We use the routine to compute a #BEMBlock of type "nearfield". */
    void CalcBlockMatrix(FlatMatrix<double> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
			 LocalHeap &lh) const override;
  };

  
  template <typename KERNEL>
  class GenericIntegralOperator : public IntegralOperator
  {
    KERNEL kernel;
  public:
    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            KERNEL _kernel,
                            struct BEMParameters param);
    
    /** Compute the sub-block of DL potential matrix which belongs to the given dofs, 
        where #matrix is dense with dim = size(testdofs) x size(trialdofs). 
        We use the routine to compute a #BEMBlock of type "nearfield". */
    void CalcBlockMatrix(FlatMatrix<double> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
			 LocalHeap &lh) const override;
  };



  template <int DIM> class LaplaceSLKernel;

  template<>
  class LaplaceSLKernel<3> 
  {
  public:
    auto Evaluate (Vec<3> x, Vec<3> y, Vec<3> nx, Vec<3> ny) const
    {
      double norm = L2Norm(x-y);
      return Mat<1,1> (1.0 / (4 * M_PI * norm));
    }
  };


  template <int DIM> class LaplaceDLKernel;

  template<>
  class LaplaceDLKernel<3> 
  {
  public:
    auto Evaluate (Vec<3> x, Vec<3> y, Vec<3> nx, Vec<3> ny) const
    {
      double norm = L2Norm(x-y);
      double nxy = InnerProduct(ny, (x-y));
      return Mat<1,1> (nxy / (4 * M_PI * norm*norm*norm));
    }
  };

}



#endif

