#ifndef NGBEM_hpp
#define NGBEM_hpp

#include "hmat.hpp"
#include "diffops.hpp"

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


  template <typename T = double>
  class IntegralOperator
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
    
    shared_ptr<HMatrix<T>> hmatrix;


  public:
    IntegralOperator (shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space, BEMParameters param);
    virtual ~IntegralOperator() = default;
    shared_ptr<BaseMatrix> GetMatrix() const { return hmatrix; }

    /** Given a #HMatrix structure, compute all block. Memory for farfield blocks gets clear here. */
    void CalcHMatrix(HMatrix<T> & hmatrix, LocalHeap &lh, struct BEMParameters &param) const;

    virtual void CalcBlockMatrix(FlatMatrix<T> matrix,
                                 FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
                                 LocalHeap &lh) const = 0;
    
    /** Compute the sub-block of integral operator matrix which belongs to the given dofs,  
        where #matrix is a #LowRankMatrix. We use the routine to compute a #BEMBlock of type "farfield". */
    virtual unique_ptr<LowRankMatrix<T>>
    CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                      LocalHeap &lh) const;
  };


  
  template <typename KERNEL>
  class GenericIntegralOperator : public IntegralOperator<typename KERNEL::value_type>
  {
    KERNEL kernel;
    typedef typename KERNEL::value_type value_type;
    typedef IntegralOperator<typename KERNEL::value_type> BASE;
    
    using BASE::trial_space; 
    using BASE::test_space;
       
    using BASE::param;

    using BASE::mapglob2bnd;
    using BASE::mapbnd2glob;
    using BASE::mapglob2bnd2;
    using BASE::mapbnd2glob2;

    using BASE::elems4dof; // contains list of elems contributing to bnd-dof
    using BASE::elems4dof2; // contains list of elems contributing to bnd-dof
    
    using BASE::trial_ct; 
    using BASE::test_ct;
    
    using BASE::hmatrix;

    shared_ptr<DifferentialOperator> trial_evaluator;
    shared_ptr<DifferentialOperator> test_evaluator;
    
  public:
    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            KERNEL _kernel,
                            struct BEMParameters param);

    GenericIntegralOperator(shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space,
                            shared_ptr<DifferentialOperator> _trial_evaluator, 
                            shared_ptr<DifferentialOperator> _test_evaluator, 
                            KERNEL _kernel,
                            struct BEMParameters param);

    
    /** Compute the sub-block of DL potential matrix which belongs to the given dofs, 
        where #matrix is dense with dim = size(testdofs) x size(trialdofs). 
        We use the routine to compute a #BEMBlock of type "nearfield". */
    void CalcBlockMatrix(FlatMatrix<value_type> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
			 LocalHeap &lh) const override;
    
    unique_ptr<LowRankMatrix<value_type>>
    CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                      LocalHeap &lh) const override;
  };



  template <int DIM> class LaplaceSLKernel;

  template<>
  class LaplaceSLKernel<3> 
  {
  public:
    typedef double value_type;
    static string Name() { return "LaplaceSL"; }

    template <typename T>        
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      return Mat<1,1,T> (1.0 / (4 * M_PI * norm));
    }
  };


  template <int DIM> class LaplaceDLKernel;

  template<>
  class LaplaceDLKernel<3> 
  {
  public:
    typedef double value_type;
    static string Name() { return "LaplaceDL"; }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      return Mat<1,1,T> (nxy / (4 * M_PI * norm*norm*norm));
    }
  };



  template <int DIM> class HelmholtzSLKernel;

  template<>
  class HelmholtzSLKernel<3> 
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzSL"; }
    
    HelmholtzSLKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      return Mat<1,1,decltype(kern)> (kern);
    }
  };


  template <int DIM> class HelmholtzDLKernel;

  template<>
  class HelmholtzDLKernel<3> 
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzDL"; }
    
    HelmholtzDLKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm);
      return Mat<1,1,decltype(kern)> (kern);
    }
  };



  template <int DIM> class CombinedFieldKernel;

  template<>
  class CombinedFieldKernel<3> 
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzDL"; }
    
    CombinedFieldKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * ( nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm)  - Complex(0,kappa)*norm*norm);
      return Mat<1,1,decltype(kern)> (kern);
    }
  };


  
}



#endif

