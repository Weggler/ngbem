#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR

#ifdef USE_FMM3D
extern "C" {
  void lfmm3d_t_c_p_(double *eps, int *nsource, double *source, double *charge,
                     int *ntarget, double *target, double *pot, int *ier);
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
    Array<Vec<3>> xpts, ypts;
  public:
    FMM_Operator(KERNEL _kernel, Array<Vec<3>> _xpts, Array<Vec<3>> _ypts)
      : kernel(_kernel), xpts(std::move(_xpts)), ypts(std::move(_ypts)) { }

    void Mult(const BaseVector & x, BaseVector & y) const override
    {
      auto fx = x.FV<typename KERNEL::value_type>();
      auto fy = y.FV<typename KERNEL::value_type>();
      
      fy = 0;
      if constexpr (std::is_same<KERNEL, class LaplaceSLKernel<3>>())
        {
#ifdef USE_FMM3D
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
#else
        for (size_t ix = 0; ix < xpts.Size(); ix++)
          for (size_t iy = 0; iy < ypts.Size(); iy++)
            fy(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);
#endif // USE_FMM3D
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

