#ifndef NGBEM_hpp
#define NGBEM_hpp



namespace ngbem
{
  using namespace ngcomp;


  class SingleLayerPotential : public SpecialElement
  {
    shared_ptr<FESpace> space;
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;

  public:
    SingleLayerPotential(shared_ptr<FESpace> aspace);

    void CalcElementMatrix(FlatMatrix<double> matrix, // matrix dim = ndof_bnd x ndof_bnd
                           LocalHeap &lh) const override;

    void GetDofNrs(Array<int> &dnums) const override;
  };




  class DoubleLayerPotential : public SpecialElement
  {
    shared_ptr<FESpace> space;
    shared_ptr<FESpace> space2;
    
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;
    Array<DofId> mapglob2bnd2;
    Array<DofId> mapbnd2glob2;

  public:
    DoubleLayerPotential(shared_ptr<FESpace> aspace, shared_ptr<FESpace> bspace);

    void CalcElementMatrix(FlatMatrix<double> matrix, // matrix dim = ndof_bnd x ndof_bnd2
                           LocalHeap &lh) const override;

    void GetDofNrs(Array<int> &dnums) const override; 
    void GetDofNrs2(Array<int> &dnums) const override;
  };

}



#endif

