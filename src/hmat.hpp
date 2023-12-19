#ifndef HMAT_hpp
#define HMAT_hpp

namespace ngbem
{
  using namespace ngcomp;

  /** Cluster. 
      A cluster is defined by its hiearchical information, 
      especially its size of dof indices and (#Number) the
      starting position in the cluster dof ordering (#PermuPos).

      A cluster is a node in the #ClusterTree. 
*/
  struct Cluster
  {
    /** Hierarchical information */
    long Level;
    long Parent, Child1, Child2;
    long Number;
    long PermuPos;
    long BasisElement;

    /** Geometrical Information */
    double Radius;
    double EVec[9], EVal[3];
    //double XMin[3],XMax[3];
    Vec<3> XMin;
    Vec<3> XMax;
    // double Centre[3];
    Vec<3> Centre;
    double DiagLength;
  };

  /** Cluster tree.
      A cluster tree is a binary tree of the index set of the FE space.
      The clustering process permutes the dofs based on their geometrical 
      distance. The geometrical information associated to a specific dof 
      may depend on the type of dof in the ho mesh.

      Each node in the tree is a #Cluster.
  */
  class ClusterTree
  {
    /** Finite element space containing the mesh. */
    shared_ptr<FESpace> space;

  public:
    /** The clustering process permutes dof indices. Therefore, we need
	a mapping from cluster indices back to FE-space dof indices. */
    Array<DofId> mapcluster2glob;

    /** Number of clusters in the tree. */
    int n_cluster;

    /** Array of length #n_cluster containing the clusters. */
    Array<Cluster> arr_clusters;

    /** Constructor (leafsize defines the minimal size of a #Cluster) */
    ClusterTree(shared_ptr<FESpace> space, int leafsize);
  };
  
  /** BEMBlock.
      A BEMBlock is a specific submatrix in a layer potential operator matrix
      in case the latter is a #HMatrix. 
      A BEMBlock can either be dense (near field) or low-rank (far field)
    
      A #HMatrix holds a list of BEMBlocks.
  */
  class BEMBlock
  {
    bool isNearField;
    Array<DofId> trialdofs;
    Array<DofId> testdofs;

    unique_ptr<BaseMatrix> matrix;
	
  public:
    /** Constructors */
    BEMBlock() {;}
    BEMBlock(FlatArray<DofId> trialdofs, FlatArray<DofId> testdofs, bool isNearField);

    /** Block type, i.e., if true it is dense and otherwise low-rank */
    bool IsNearField() const {return isNearField; }
    
    /** Access to matrix */  
    const unique_ptr<BaseMatrix> & GetMat() const { return matrix; }
    void SetMat(unique_ptr<BaseMatrix> _matrix) { matrix = std::move(_matrix); }
    
    /** All dofs of the trial functions contributing to the block */  
    FlatArray<DofId> GetTrialDofs() const { return trialdofs; }
    /** All dofs of the test functions contributing to the block */  
    FlatArray<DofId> GetTestDofs() const { return testdofs; }

    /** Matrix-vector-multiply with global vectors x and y. Only the entries specified
	by #testdofs and #trialdofs are used. */
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const;
    void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;    
  };

  /** Low-rank matrix
      A low-rank matrix \f$C\f$ of dimension \$fn\times m\$f has the form
      \f$ C = A B^\top \f$, where \f$A, B\f$ are matrices of dimensions
      \$fn\times r\$f and \f$m \times r \f$ with \$fr<m,n\f$.
  */

  template <typename T = double>
  class LowRankMatrix : public BaseMatrix 
  {
    /** Height of the matrix #A. */
    // size_t m;   // directly ask A
    /** Width of the matrix #Bt. */
    // size_t n;
    /** The rank is the width of #A and height of #Bt. */
    size_t rank;
    /** The Matrix \f$A\f$ of the low-rank factorization. */
    Matrix<T> A;
    /** The transpose of the matrix \f$B\f$ of the low-rank factorization. */
    Matrix<T> Bt;
  public:
    /** Empty constructor. */
    // LowRankMatrix();
    /** Constructor. */
    LowRankMatrix(Matrix<T> A, Matrix<T> Bt);

    /** Access matrices of factorization . */
    const auto & GetA() const { return A; }
    const auto & GetBt() const { return Bt; }
    
    /** Matrix-vector-multiplication, e.g. \f$y = y + s (A B^\top ) y \f$. */
    virtual void Mult(const BaseVector & x, BaseVector & y) const override;    
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultAdd (Complex s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const override;

    bool IsComplex() const override { return false; }

    virtual int VHeight() const override { return A.Height(); }
    virtual int VWidth() const override { return Bt.Width(); }
    virtual int VRank() const { return rank; }

    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
  };

  /** Hierarchical matrix.
      A hierarchical matrix is a generalized block matrix with a partition
      based on clusters. It consists of near-field dense blocks and far-field
      low-rank blocks. A block is  #BEMBlock. 
  */

  template <typename T = double>
  class HMatrix : public BaseMatrix
  {
    /** Row cluster tree. */
    shared_ptr<ClusterTree> row_ct;
    
    /** Column cluster tree. */
    shared_ptr<ClusterTree> col_ct;
    
    /** Parameter in admissibility condition. Decides whether a block is part
	of the near- or far-field. */
    const double eta;
    
    /** The hierarchical matrix realized as list of admissible and inadmissible blocks. */
    Array<BEMBlock> matList;

    /** HACK: the solver routines need the number of all dofs, not only boundary dofs. */
    int width_vol_dof, height_vol_dof;
    
  public:
    HMatrix(shared_ptr<ClusterTree> col_ct, shared_ptr<ClusterTree> row_ct, double eta, int width_vol_dof, int height_vol_dof);

    Array<BEMBlock> & GetMatList() { return matList; }

    /** Matrix-vector-multiply-add: \f$y = y + s A x \f$. */
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;
    bool IsComplex() const override { return false; }

    /*
    virtual int VHeight() const override { return row_ct->mapcluster2glob.Size(); }
    virtual int VWidth() const override { return col_ct->mapcluster2glob.Size(); }
    */
    virtual int VHeight() const override { return height_vol_dof; }
    virtual int VWidth() const override { return width_vol_dof; }
    
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;

    /** Get the size of the allocated storage of the hierarchical matrix in bytes. */
    //size_t GetMemSize() const;
  };
}

#endif

