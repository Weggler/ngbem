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




  // rotated 
  class DiffOpMaxwell : public DiffOp<DiffOpMaxwell>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 3 };
    enum { DIM_ELEMENT = 2 };
    enum { DIM_DMAT = 4 };
    enum { DIFFORDER = 1 };

    static string Name() { return "Maxwell"; }
    
    static const HCurlFiniteElement<2> & Cast (const FiniteElement & fel) 
    { return static_cast<const HCurlFiniteElement<2>&> (fel); }

    ///
    // mat is 4 x ndof
    template <typename AFEL, typename MIP, typename MAT>
    static void GenerateMatrix (const AFEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    {
      // mat.AddSize(4, fel.GetNDof()) = 0.0;
      auto matvec = mat.Rows(0,3);
      Cast(fel).CalcMappedShape (mip, Trans(matvec));
      for (int i = 0; i < fel.GetNDof(); i++)
        {
          Vec<3> shape = matvec.Col(i);
          matvec.Col(i) = Cross(mip.GetNV(), shape);
        }
      // *testout << "mat1 = " << mat << endl;
      // mat.AddSize(4, fel.GetNDof()) = 0.0;      
      mat.Row(3) =
        1.0/mip.GetJacobiDet() * 
	Cast(fel).GetCurlShape(mip.IP(),lh).Col(0);
      // *testout << "scalar mat = " << endl << mat << mat;
      // *testout << "mat2 = " << mat << endl;      
    }

    /// mat is (ndof*4) x mip.Size()
    static void GenerateMatrixSIMDIR (const FiniteElement & fel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      Cast(fel).CalcMappedShape (mir, mat.Rows(0, 3*fel.GetNDof()));
      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<2,3>>&>(mir[j]).GetNV();
          for (int i = fel.GetNDof()-1; i >= 0; i--)
            {
              Vec<3,SIMD<double>> shape = mat.Col(j).Range(3*i, 3*i+3);
              mat.Col(j).Range(4*i,4*i+3) = Cross(nv, shape);
              // mat.Col(j).Range(4*i+3,4*i+4) = SIMD<double>(0.0);
            }
        }
      // *testout << "simd mat1 = " << endl << mat.AddSize(4*fel.GetNDof(), mir.Size()) << endl;            
      // mat.AddSize(4*fel.GetNDof(), mir.Size()) = SIMD<double>(0.0);
      
      constexpr size_t BS=16;
      LocalHeapMem<BS*SIMD<double>::Size()*sizeof(SIMD<MappedIntegrationPoint<2,2>>)+64> lh("genmatlh");
      FE_ElementTransformation<2,2> trafo2d(fel.ElementType());
      for (size_t first = 0; first < mir.Size(); first += BS)
        {
          HeapReset hr(lh);
          size_t next = std::min(first+BS, mir.Size());
          SIMD_MappedIntegrationRule<2,2> mir2d(mir.IR().Range(first, next), trafo2d, lh);
          Cast(fel).CalcMappedCurlShape (mir2d, mat.RowSlice(3,4).Cols(first, next));
        }
      for (size_t i = 0; i < mir.Size(); i++)
        mat.Col(i).Slice(3,4).Range(fel.GetNDof()) *= 1.0 / mir[i].GetJacobiDet();
      // *testout << "simd mat2 = " << endl << mat.AddSize(4*fel.GetNDof(), mir.Size()) << endl;      
    }

    using DiffOp<DiffOpMaxwell>::ApplySIMDIR;
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<Complex> x, BareSliceMatrix<SIMD<Complex>> y)
    {
      Cast(fel).Evaluate (mir, x, y.Rows(3));
      for (int j = 0; j < mir.Size(); j++)
        {
          Vec<3,SIMD<double>> nv = static_cast<const SIMD<ngfem::MappedIntegrationPoint<2,3>>&>(mir[j]).GetNV();
          Vec<3,SIMD<Complex>> shape = y.Col(j).Range(0,3);
          y.Col(j).Range(0,3) = Cross(nv, shape);
        }
      y.Row(3).Range(mir.Size()) = SIMD<Complex>(0.0);
      // TODO: curl part
    }
  };
  
}


#endif
