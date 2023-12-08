#ifndef HMAT_hpp
#define HMAT_hpp

namespace ngbem
{
  using namespace ngcomp;

  /** Cluster. */
  struct Cluster
  {
    long Level;
    long Parent, Child1, Child2;
    long Number;
    long PermuPos;
    long BasisElement;

    double Radius;
    double EVal[3],EVec[9];
    double XMin[3],XMax[3];
    double Centre[3];
    double DiagLength;
  };

  /** Cluster tree.
      A cluster tree is a binary tree of the index set of the FE space associated
      with its geometrical information, i.e. the physical support points of the
      finite element functions on the mesh.

      A cluster in the tree is defined by its size (number of indices) and its
      starting position in the cluster ordering.
  */
  class ClusterTree
  {
    /** Finite element spaces containing the mesh. */
    shared_ptr<FESpace> space;

  public:
    /** The clustering process only consider boundary dofs. Therefore, we need
	a mapping from cluster indices to global indices. */
    Array<DofId> mapcluster2glob;

    /** Number of clusters in the tree. */
    int n_cluster;

    /** Array of length #n_cluster containing the clusters. */
    Array<struct Cluster> arr_clusters;

    ClusterTree(shared_ptr<FESpace> space, int leafsize);
  };

  class BEMBlock
  {
    bool isNearField;
    Array<DofId> setI;
    Array<DofId> setJ;

    shared_ptr<BaseMatrix> matrix;
	
  public:
    BEMBlock() {;}
    BEMBlock(Array<DofId> &setI, Array<DofId> &setJ, bool isNearField);
    bool IsNearField() const {return isNearField; }
    shared_ptr<shared_ptr<BaseMatrix>> GetMat() const { return make_shared<shared_ptr<BaseMatrix>>(matrix); }
    shared_ptr<Array<DofId>> GetI() const { return make_shared<Array<DofId>>(setI); }
    shared_ptr<Array<DofId>> GetJ() const { return make_shared<Array<DofId>>(setJ); }

    /** Matrix-vector-multiply with global vectors x and y. Only the entries specified
	by #setI and #setJ are used. */
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;    
  };

  /** Low-rank matrix
      A low-rank matrix \f$C\f$ of dimension \$fn\times m\$f has the form
      \f$ C = A B^\top \f$, where \f$A, B\f$ are matrices of dimensions
      \$fn\times r\$f and \f$m \times r \f$ with \$fr<m,n\f$.
  */
  class LowRankMatrix : public BaseMatrix 
  {
    /** Height of the matrix #A. */
    size_t m;
    /** Width of the matrix #Bt. */
    size_t n;
    /** The rank is the width of #A and height of #Bt. */
    size_t rank;
    /** The Matrix \f$A\f$ of the low-rank factorization. */
    shared_ptr<Matrix<>> A;
    /** The transpose of the matrix \f$B\f$ of the low-rank factorization. */
    shared_ptr<Matrix<>> Bt;
  public:
    /** Empty constructor. */
    LowRankMatrix();
    /** Full cunstructor. */
    LowRankMatrix(shared_ptr<Matrix<>> A, shared_ptr<Matrix<>> Bt);
    
    /** Matrix-vector-multiply-add: \f$y = y + s A B^\top y \f$. */
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    shared_ptr<Matrix<>> GetA() const { return A; }
    shared_ptr<Matrix<>> GetBt() const { return Bt; }

    
    bool IsComplex() const  { return false; }

    virtual int VHeight() const { return m; }
    virtual int VWidth() const { return n; }
    virtual int VRank() const { return rank; }

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
  };

  /** Hierarchical matrix.
      A hierarchical matrix is a generalized block matrix with a partition
      based on clusters. It consists of near-field dense blocks and far-field
      low-rank blocks.
  */
  class HMatrix : public BaseMatrix
  {
    /** Row cluster tree. */
    shared_ptr<ClusterTree> row_ct;
    
    /** Column cluster tree. */
    shared_ptr<ClusterTree> col_ct;
    
    /** Parameter in admissibility condition. Decides whether a block is part
	of the near- or far-field.
    */
    const double eta;
    
    /** The hierarchical matrix realized as list of admissible and inadmissible blocks. */
    Array<BEMBlock> matList;

    /** HACK: the solver routines need the number of all dofs, not only boundary dofs. */
    int width_vol_dof, height_vol_dof;
    
  public:
    HMatrix(shared_ptr<ClusterTree> row_ct, shared_ptr<ClusterTree> col_ct, double eta, int width_vol_dof, int height_vol_dof);

    shared_ptr<Array<BEMBlock>> GetMatList() const { return make_shared<Array<BEMBlock>>(matList); }

    /** Matrix-vector-multiply-add: \f$y = y + s A B^\top x \f$. */
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    bool IsComplex() const { return false; }

    virtual int VHeight() const { return row_ct->mapcluster2glob.Size(); }
    virtual int VWidth() const { return col_ct->mapcluster2glob.Size(); }

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
  };
}

#endif

