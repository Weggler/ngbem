#ifndef FILE_FMMOPERATOR
#define FILE_FMMOPERATOR



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
      auto fx = x.FV<double>();
      auto fy = y.FV<double>();

      fy = 0;
      if (std::is_same<KERNEL, class LaplaceSLKernel<3>>())
        for (size_t ix = 0; ix < xpts.Size(); ix++)
          for (size_t iy = 0; iy < ypts.Size(); iy++)
            fy(iy) += __Kernel(xpts[ix], ypts[iy]) * fx(ix);
      else
        throw Exception("fmm not available");
      // fmmapply
    }


    AutoVector CreateRowVector () const override
    {
      return make_unique<VVector<double>>(xpts.Size());
    }
    AutoVector CreateColVector () const override
    {
      return make_unique<VVector<double>>(ypts.Size());
    }
  };
}

#endif

