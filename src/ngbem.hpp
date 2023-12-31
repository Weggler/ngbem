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


  /** The IntegralOperator provides methods for the assembly of and access to a matrix 
      resulting from a variational formulation of a boundary integral equation.*/
  template <typename T = double>
  class IntegralOperator
  {
  protected:
    shared_ptr<FESpace> trial_space; 
    shared_ptr<FESpace> test_space; 
    
    /* parameters that specify the appoximation of linear operators and solving */
    BEMParameters param;

    /* boundary to global fe dofs mappings */
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;
    Array<DofId> mapglob2bnd2;
    Array<DofId> mapbnd2glob2;

    Table<int> elems4dof; // contains list of elems contributing to bnd-dof
    Table<int> elems4dof2; // contains list of elems contributing to bnd-dof
    
    /* #ClusterTree define cluster pairs, i.e., the blocks of the hmatrix */
    shared_ptr<ClusterTree> trial_ct; 
    shared_ptr<ClusterTree> test_ct;
    
    shared_ptr<HMatrix<T>> hmatrix;


  public:

    /** Constructor. */
    IntegralOperator (shared_ptr<FESpace> _trial_space, shared_ptr<FESpace> _test_space, BEMParameters param);
    virtual ~IntegralOperator() = default;

    /** GetHMatrix returns the #hmatrix. */
    shared_ptr<BaseMatrix> GetMatrix() const { return hmatrix; }

    /** CalcHMatrix fills the #hmatrix, i.e., memory for farfield blocks is allocated and 
        all blocks are computed. */
    void CalcHMatrix(HMatrix<T> & hmatrix, LocalHeap &lh, struct BEMParameters &param) const;

    /** CalcBlockMatrix computes the block of entries with trialdofs and testdofs indices. */
    virtual void CalcBlockMatrix(FlatMatrix<T> matrix,
                                 FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
                                 LocalHeap &lh) const = 0;
    
    /** CalcFarFieldBlock computes a low-rank approximation of block with trialdofs and testdofs. */
    virtual unique_ptr<LowRankMatrix<T>>
    CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                      LocalHeap &lh) const;
  };


  
  /** The GenericIntegralOperator is a templated #IntegralOperator, the template type is 
      the kernel the specific potential,i.e. a fundamental solution or its 
      derivative of specific pde. */
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

    using BASE::elems4dof; 
    using BASE::elems4dof2; 
    
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

    void CalcBlockMatrix(FlatMatrix<value_type> matrix, FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, 
			 LocalHeap &lh) const override;
    
    unique_ptr<LowRankMatrix<value_type>>
    CalcFarFieldBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs,
                      LocalHeap &lh) const override;
  };




  /** LaplaceSLkernel is the kernel for the single layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM> class LaplaceSLKernel;

  /** LaplaceSLkernel in 3D reads 
      $$ G(x-y) = \frac{1}{4\,\pi \, | x-y| }, \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
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
      return 1.0 / (4 * M_PI * norm);   
      // return Mat<1,1,T> (1.0 / (4 * M_PI * norm));
    }
  };


  /** LaplaceDLkernel is the kernel for the double layer potential of 
      the Laplace equation $ \Delta u = 0 \,.$  */
  template <int DIM> class LaplaceDLKernel;

  /** LaplaceDLkernel in 3D reads 
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, 
          \frac{ \langle n(y), x-y\rangle }{ | x-y|^3 }, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
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
      return nxy / (4 * M_PI * norm*norm*norm);
      // return Mat<1,1,T> (nxy / (4 * M_PI * norm*norm*norm));
    }
  };



  /** HelmholtzSLkernel is the kernel for the double layer potential of the 
      Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,. $ */
  template <int DIM> class HelmholtzSLKernel;

  /** HelmholtzSLkernel in 3D reads 
      $$ G(x-y) = \frac{1 }{4\,\pi} \,\frac{e^{i\,\kappa \, |x-y| }{|x-y|} \, 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class HelmholtzSLKernel<3> 
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "HelmholtzSL"; }
    
    /** Construction of the kernel specifies the wavenumber $\kappa$. */
    HelmholtzSLKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm);
      return kern;
      // return Mat<1,1,decltype(kern)> (kern);
    }
  };


  /** HelmholtzSLkernel is the kernel for the double layer potential of 
      the Helmholtz equation $ -\Delta u - \kappa^2 u = 0, \; \kappa>0\,.$ */
  template <int DIM> class HelmholtzDLKernel;

  /** HelmholtzSLkernel in 3D reads
      $$ \frac{\partial }{ \partial n_y} G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \left( 1 - i\,\kappa\, | x-y| \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
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
      return kern;
      // return Mat<1,1,decltype(kern)> (kern);
    }
  };



  /** CombinedFieldKernel is a kernel for the combined field integral equation 
      is considered for the Helmholtz equation. */
  template <int DIM> class CombinedFieldKernel;

  /** CombinedFildKernel in 3D reads
      $$ G(x-y) = \frac{1}{4\,\pi} \, \frac{e^{i\,\kappa\,|x-y|}}{|x-y|^3} \, 
          \left( \langle n_y, x-y\rangle - i\,\kappa\, | x-y| - i\,\kappa\,|x-y|^2 \right), 
          \quad x, y \in \mathbb R^3, \; x\not=y\,. $$ */
  template<>
  class CombinedFieldKernel<3> 
  {
    double kappa;
  public:
    typedef Complex value_type;
    static string Name() { return "Helmholtz Combined Field"; }
    
    CombinedFieldKernel (double _kappa) : kappa(_kappa) { }

    template <typename T>    
    auto Evaluate (Vec<3,T> x, Vec<3,T> y, Vec<3,T> nx, Vec<3,T> ny) const
    {
      T norm = L2Norm(x-y);
      T nxy = InnerProduct(ny, (x-y));
      auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
        * ( nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm)  - Complex(0,kappa)*norm*norm);
      return kern;
      // return Mat<1,1,decltype(kern)> (kern);
    }
  };


  
}



#endif

