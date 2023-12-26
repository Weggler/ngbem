#ifndef DIFFOPS_HPP
#define DIFFOPS_HPP

#include <bla.hpp>

namespace ngbem
{
  using namespace ngbla;



  class DiffOpBoundaryRot : public DiffOp<DiffOpBoundaryRot>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 1 };

    static string Name() { return "boundaryrot"; }
    
    static const ScalarFiniteElement<2> & Cast (const FiniteElement & fel) 
    { return static_cast<const ScalarFiniteElement<2>&> (fel); }

    ///
    // mat is 3 x ndof
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedDShape (mip, Trans(mat));
      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Vec<3> grad = mat.Col(i);
          mat.Col(i) = Cross(mip.GetNV(), grad);
        }
    }

    /// mat is (ndof*3) x mip.Size()
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedDShape (mir, mat);

      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<3,3>>&>(mir[j]).GetNV();
          for (int i = 0; i < fel.GetNDof(); i++)
            {
              Vec<3,SIMD<double>> grad = mat.Col(j).Range(3*i, 3*i+3);
              mat.Col(j).Range(3*i,3*i+3) = Cross(nv, grad);
            }
        }
    }
  };




  class DiffOpRotatedTrace : public DiffOp<DiffOpRotatedTrace>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 3 };
    enum { DIFFORDER = 1 };

    static string Name() { return "rotatedtrace"; }
    
    static const HCurlFiniteElement<2> & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlFiniteElement<2>&> (fel); }

    ///
    // mat is 3 x ndof
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      Cast(fel).CalcMappedShape (mip, Trans(mat));
      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Vec<3> shape = mat.Col(i);
          mat.Col(i) = Cross(mip.GetNV(), shape);
        }
    }

    /// mat is (ndof*3) x mip.Size()
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedShape (mir, mat);

      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<3,3>>&>(mir[j]).GetNV();
          for (int i = 0; i < fel.GetNDof(); i++)
            {
              Vec<3,SIMD<double>> shape = mat.Col(j).Range(3*i, 3*i+3);
              mat.Col(j).Range(3*i,3*i+3) = Cross(nv, shape);
            }
        }
    }
  };
  
}


#endif
