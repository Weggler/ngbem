#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR

#ifdef USE_KiFMM
// namespace kifmm {
  #include <kifmm_rs.h>
// }
#endif // USE_KiFMM

#ifdef USE_FMM3D
extern "C" {
  void lfmm3d_t_c_p_(double *eps, int *nsource, double *source, double *charge,
                     int *ntarget, double *target, double *pot, int *ier);
  void hfmm3d_t_c_p_(double* eps, std::complex<double>* zk,
                     int* nsources, double* sources,
                     std::complex<double>* charges, int* ntargets, double* targets,
                     std::complex<double>* potentials, int* ier);
  void hfmm3d_t_cd_p_(double* eps, std::complex<double>* zk,
                      int* nsources, double* sources,
                      std::complex<double>* charges,
                      std::complex<double>* dipvec,
                      int* ntargets, double* targets,
                      std::complex<double>* potentials, int* ier);
}
#endif // USE_FMM3D

namespace ngbem
{
  
  inline double __Kernel(Vec<3> px, Vec<3> py)
  {
    if (L2Norm2(px-py)==0) return 0.0;
    return 1.0 / (4*M_PI) / L2Norm(px-py);
  }


  template <typename TSCAL>
  class Base_FMM_Operator : public BaseMatrix
  {
  protected:
    Array<Vec<3>> xpts, ypts, xnv, ynv;
  public:
    Base_FMM_Operator(Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                   Array<Vec<3>> _xnv, Array<Vec<3>> _ynv)
      : xpts(std::move(_xpts)), ypts(std::move(_ypts)),
        xnv(std::move(_xnv)), ynv(std::move(_ynv))
    { }
    
    AutoVector CreateRowVector () const override
    {
      return make_unique<VVector<TSCAL>>(xpts.Size());
    }
    AutoVector CreateColVector () const override
    {
      return make_unique<VVector<TSCAL>>(ypts.Size());
    }
  };


  
  template <typename KERNEL>
  class FMM_Operator : public Base_FMM_Operator<typename KERNEL::value_type>
  {
    KERNEL kernel;
    typedef Base_FMM_Operator<typename KERNEL::value_type> BASE;
    using BASE::xpts, BASE::ypts, BASE::xnv, BASE::ynv;

  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                 Array<Vec<3>> _xnv, Array<Vec<3>> _ynv)
      : BASE(std::move(_xpts), std::move( _ypts), std::move(_xnv), std::move(_ynv)),
      kernel(_kernel)
    { }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply"); RegionTimer reg(tall);
      
      auto fx = x.FV<typename KERNEL::value_type>();
      auto fy = y.FV<typename KERNEL::value_type>();
      
      fy = 0;
      if constexpr (std::is_same<KERNEL, class LaplaceSLKernel<3>>())
        {
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              fy(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);
        }
      else if constexpr (std::is_same<KERNEL, class HelmholtzSLKernel<3>>())
        {
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              {
                double norm = L2Norm(xpts[ix]-ypts[iy]);
                if (norm > 0)
                  {
                    auto kern = exp(Complex(0,kernel.GetKappa())*norm) / (4 * M_PI * norm);
                    fy(iy) += kern * fx(ix);
                  }
              }
        }
      else if constexpr (std::is_same<KERNEL, class CombinedFieldKernel<3>>())
        {
          /*
            T norm = L2Norm(x-y);
            T nxy = InnerProduct(ny, (x-y));
            auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
            * ( nxy * (Complex(1,0)*T(1.) - Complex(0,kappa)*norm)  - Complex(0,kappa)*norm*norm);
            // return kern;
            return Vec<1,decltype(kern)> (kern);
          */
          double kappa = kernel.GetKappa();
          for (size_t ix = 0; ix < xpts.Size(); ix++)
            for (size_t iy = 0; iy < ypts.Size(); iy++)
              {
                double norm = L2Norm(xpts[ix]-ypts[iy]);
                if (norm > 0)
                  {
                    // auto kern = exp(Complex(0,kernel.GetKappa())*norm) / (4 * M_PI * norm);
                    // fy(iy) += kern * fx(ix);
                    double nxy = InnerProduct(ynv[iy], xpts[ix]-ypts[iy]);
                    auto kern = exp(Complex(0,kappa)*norm) / (4 * M_PI * norm*norm*norm)
                      * ( nxy * (1.0 - Complex(0,kappa)*norm) - Complex(0,kappa)*norm*norm);
                    fy(iy) += kern * fx(ix);                    
                  }
              }
        }
      else
        throw Exception("fmm not available");
    }

  };

  

  
#ifdef USE_FMM3D
  
  // we can specialze the complete class ... 
  template <> 
  class FMM_Operator<LaplaceSLKernel<3>> : public Base_FMM_Operator<double>
  {
    typedef LaplaceSLKernel<3> KERNEL;
    typedef Base_FMM_Operator<double> BASE;
    
    KERNEL kernel;
    using BASE::xpts, BASE::ypts, BASE::xnv, BASE::ynv;    

  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                 Array<Vec<3>> _xnv, Array<Vec<3>> _ynv)
      : BASE(std::move(_xpts), std::move( _ypts), std::move(_xnv), std::move(_ynv)),
        kernel(_kernel)
    { }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply LaplaceSL (FMM3D)"); RegionTimer reg(tall);
      auto fx = x.FV<double>();
      auto fy = y.FV<double>();
      
      fy = 0;

      double eps = 1e-8;
      int ier;
      int size_x = xpts.Size();
      int size_y = ypts.Size();
      lfmm3d_t_c_p_(&eps, &size_x, xpts[0].Data(),
                    fx.Data(), &size_y, ypts[0].Data(),
                    fy.Data(), &ier);
      if (ier != 0)
        throw Exception("FMM3D failed with err code " + std::to_string(ier));
      y *= 1.0 / (4*M_PI);
    }
  };

  
  // or just specialze the apply method:

  template <>
  void FMM_Operator<HelmholtzSLKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply HelmholtzSL (FMM3D)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();

    fy = 0;

    double eps = 1e-8;
    int ier;
    std::complex<double> zk = kernel.GetKappa();
    int size_x = xpts.Size();
    int size_y = ypts.Size();
    hfmm3d_t_c_p_(&eps, &zk, &size_x, xpts[0].Data(),
                  fx.Data(), &size_y, ypts[0].Data(),
                  fy.Data(), &ier);
    if (ier != 0)
      throw Exception("FMM3D failed with err code " + std::to_string(ier));
    y *= 1.0 / (4*M_PI);
  }

  
  template <>
  void FMM_Operator<CombinedFieldKernel<3>> :: Mult(const BaseVector & x, BaseVector & y) const 
  {
    static Timer tall("ngbem fmm apply CombinedField (FMM3D)"); RegionTimer reg(tall);
    auto fx = x.FV<Complex>();
    auto fy = y.FV<Complex>();
    
    fy = 0;
    
    double eps = 1e-8;
    int ier;
    Complex zk = kernel.GetKappa();
    Vector<Complex> c(xpts.Size());
    for(auto i : Range(c))
      c[i] = -Complex(0, kernel.GetKappa()) * fx[i];
    Vector<Vec<3, Complex>> v(xpts.Size());
    for(auto i : Range(xnv))
      v[i] = fx[i] * xnv[i];
    int size_x = xpts.Size();
    int size_y = ypts.Size();
    hfmm3d_t_cd_p_(&eps, &zk, &size_x, xpts[0].Data(),
                   c.Data(), v[0].Data(),
                   &size_y, ypts[0].Data(),
                   fy.Data(), &ier);
    if (ier != 0)
      throw Exception("FMM3D failed with err code " + std::to_string(ier));
    for(auto i : Range(fy))
      fy[i] *= 1 / (4*M_PI);
  }

#endif // USE_FMM3D


#ifdef USE_KiFMM

  // we can specialze the complete class ... 
  template <> 
  class FMM_Operator<LaplaceSLKernel<3>> : public Base_FMM_Operator<double>
  {
    typedef LaplaceSLKernel<3> KERNEL;
    typedef Base_FMM_Operator<double> BASE;
    
    KERNEL kernel;
    using BASE::xpts, BASE::ypts, BASE::xnv, BASE::ynv;

    // LaplaceFft64* fmm;
    // Array<size_t> expansion_order;

  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                 Array<Vec<3>> _xnv, Array<Vec<3>> _ynv)
      : BASE(std::move(_xpts), std::move( _ypts), std::move(_xnv), std::move(_ynv)),
        kernel(_kernel)
    {
      // do some initialization ? 
    }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply LaplaceSL (FMM3D)"); RegionTimer reg(tall);
      auto fx = x.FV<double>();
      auto fy = y.FV<double>();

      // following: https://github.com/bempp/kifmm/blob/enh/c-abi/kifmm/c/main.c
      bool prune_empty = true;
      uint64_t n_crit = 150;
      uint64_t depth = 0;
      double singular_value_threshold = 0.001;
      
      uintptr_t expansion_order[] = {6};
      uintptr_t nexpansion_order = 1;
      uintptr_t block_size = 32;

      // Instantiate a Laplace evaluator
      struct FmmEvaluator *evaluator =
        laplace_fft_f64_alloc(
                              expansion_order, nexpansion_order,
                              xpts[0].Data(), 3*xpts.Size(), 
                              ypts[0].Data(), 3*ypts.Size(), 
                              fx.Data(), fx.Size(), 
                              prune_empty, n_crit, depth, block_size
                              );
      
      bool timed = true;
      evaluate(evaluator, timed);

      MortonKeys *leaves = leaves_target_tree(evaluator);
      
      // printf("Number of leaf keys (n): %zu\n", leaves->len);
      // for (uintptr_t i = 0; i < 5; ++i) {
      // printf("Element %zu: %llu\n", i, leaves->data[i]);
      // }

      Potentials *potentials = leaf_potentials(evaluator, leaves->data[123]);

      /*
      Coordinates *coordinates =
        coordinates_source_tree(evaluator, leaves->data[123]);
      printf("Number of coordinates: %zu\n", coordinates->len);

      const double *coords = (const double *)coordinates->data;
      for (uintptr_t i = 0; i < 5; ++i) {
        printf("Element %zu: [%f, %f, %f]\n", i, coords[i * 3], coords[i * 3 + 1],
               coords[i * 3 + 2]);
      }
      */
      
      // printf("Number of potentials: %zu\n", potentials->n);
      Potential *pot = &potentials->data[0];
      const double *data = (const double *)pot->data;
      for (uintptr_t i = 0; i < 5; ++i) {
        // printf("Element %zu: %f\n", i, data[i]);
        fy(i) = data[i];
      }

      // Cleanup
      free_fmm_evaluator(evaluator);
      free_morton_keys(leaves);
      free(potentials);


      
      /*
      fy = 0;

      Array<size_t> expansion_order =
        { 5 };
      auto fmm = laplace_blas_svd_f64_alloc
        (expansion_order.Data(), expansion_order.Size(),
         xpts[0].Data(), xpts.Size() * 3,
         ypts[0].Data(), ypts.Size() * 3,
         fx.Data(), fx.Size(), // charges
         true,                       // prune empty
         150,                        // n_crit
         0,                          // depth
         1e-8);
      evaluate(fmm, false);
      auto potentials = potentials(fmm);
      fy = FlatVector<double>(potentials.len, (double*) potentials.data);
      */
      
    }
  };

#endif // KiFMM
  
  
}


#endif

