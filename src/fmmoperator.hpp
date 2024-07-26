#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR

#ifdef USE_KiFMM
namespace kifmm {
  #include <kifmm_rs.h>
}
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

  template <typename KERNEL>
  class FMM_Operator : public BaseMatrix
  {
    KERNEL kernel;
    Array<Vec<3>> xpts, ypts, xnv, ynv;
#ifdef USE_KiFMM
    kifmm::LaplaceFft64* fmm;
    Array<size_t> expansion_order;
#endif // USE_KiFMM

  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts,
                 Array<Vec<3>> _xnv, Array<Vec<3>> _ynv)
      : kernel(_kernel), xpts(std::move(_xpts)), ypts(std::move(_ypts)),
        xnv(std::move(_xnv)), ynv(std::move(_ynv))
    {
    }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      static Timer tall("ngbem fmm apply"); RegionTimer reg(tall);
      
      auto fx = x.FV<typename KERNEL::value_type>();
      auto fy = y.FV<typename KERNEL::value_type>();
      
      fy = 0;
      if constexpr (std::is_same<KERNEL, class LaplaceSLKernel<3>>())
        {
#ifdef USE_KiFMM
          static Timer t("KiFMM Apply");
          RegionTimer r(t);
          Array<size_t> expansion_order =
            { 5 };
          auto fmm = kifmm::laplace_blas_svd_f64_alloc
            (expansion_order.Data(), expansion_order.Size(),
             xpts[0].Data(), xpts.Size() * 3,
             ypts[0].Data(), ypts.Size() * 3,
             fx.Data(), fx.Size(), // charges
             true,                       // prune empty
             150,                        // n_crit
             0,                          // depth
             1e-8);
          kifmm::evaluate(fmm, false);
          auto potentials = kifmm::potentials(fmm);
          fy = FlatVector<double>(potentials.len, (double*) potentials.data);
#elif USE_FMM3D
          static Timer t("nbem - FMM3D Apply Laplace");
          RegionTimer r(t);
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
#else // USE_FMM3D (also no KiFMM)
        for (size_t ix = 0; ix < xpts.Size(); ix++)
          for (size_t iy = 0; iy < ypts.Size(); iy++)
            fy(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);
#endif // USE_KiFMM

        }
      else if constexpr (std::is_same<KERNEL, class HelmholtzSLKernel<3>>())
        {
#ifdef USE_KiFMM
#elif USE_FMM3D
          static Timer t("nbem - FMM3D Apply Helmholtz");
          RegionTimer r(t);
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
#else // USE_FMM3D
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
#endif // USE_FMM3D
        }


      else if constexpr (std::is_same<KERNEL, class CombinedFieldKernel<3>>())
        {
#ifdef USE_KiFMM
          // todo
#elif USE_FMM3D
          static Timer t("nbem - FMM3D Apply Helmholtz");
          RegionTimer r(t);
          double eps = 1e-8;
          int ier;
          std::complex<double> zk = kernel.GetKappa();
          Vector<std::complex<double>> c(xpts.Size());
          for(auto i : Range(c))
            c[i] = -Complex(0, kernel.GetKappa()) * fx[i];
          Vector<Vec<3, std::complex<double>>> v(xpts.Size());
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
#else // USE_FMM3D
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
#endif // USE_FMM3D || USE_KiFMM
        }
      else
        throw Exception("fmm not available");
    }


    AutoVector CreateRowVector () const override
    {
      return make_unique<VVector<typename KERNEL::value_type>>(xpts.Size());
    }
    AutoVector CreateColVector () const override
    {
      return make_unique<VVector<typename KERNEL::value_type>>(ypts.Size());
    }
  };
}

#endif

