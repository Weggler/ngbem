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
    /** The clustering produces a different ordering of boundary dofs. */
    Array<DofId> mapbnd2cluster;

    /** Number of clusters in the tree. */
    int n_cluster;

    /** Array of length #n_cluster containing the clusters. */
    Array<struct Cluster> arr_clusters;

    ClusterTree(shared_ptr<FESpace> space, int leafsize);
  };

  /** Hierarchical matrix.
      A hierarchical matrix is a generalized block matrix with a partition
      based on clusters. It consists of near-field dense blocks and far-field
      low-rank blocks.
   */
  class HMatrix
  {
    /** Row cluster tree. */
    shared_ptr<ClusterTree> row_ct;
    
    /** Column cluster tree. */
    shared_ptr<ClusterTree> col_ct;
    
    /** Parameter in admissibility condition. Decides whether a block is part
	of the near- or far-field.
     */
    const double eta;
    
    /** The hierarchical matrix realized as a BlockMatrix. */
    shared_ptr<BaseMatrix> mat;
    
  public:
    HMatrix(shared_ptr<ClusterTree> row_ct, shared_ptr<ClusterTree> col_ct, double eta);

    /** Matrix-vector-multiply-add: \f$y = y + s A B^\top y \f$. */
    //void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;

    bool IsComplex() const { return false; }
    //virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override;

    virtual int VHeight() const { return row_ct->mapbnd2cluster.Size(); }
    virtual int VWidth() const { return col_ct->mapbnd2cluster.Size(); }

    //virtual AutoVector CreateRowVector () const;
    //virtual AutoVector CreateColVector () const;
  };

  /** Low-rank matrix
      A low-rank matrix \f$C\f$ of dimension \$fn\times m\$f has the form
      \f$ C = A B^\top \f$, where \f$A, B\f$ are matrices of dimensions
      \$fn\times r\$f and \f$m \times r \f$ with \$fr<m,n\f$.
   */
  class LowRankMatrix 
  {
    /** Height of the matrix #A. */
    size_t m;
    /** Height of the matrix #Bt. */
    size_t n;
    /** The rank is the width of #A and #Bt. */
    size_t rank;
    /** The Matrix \f$A\f$ of the low-rank factorization. */
    shared_ptr<Matrix<>> A;
    /** The transpose of the matrix \f$B\f$ of the low-rank factorization. */
    shared_ptr<Matrix<>> Bt;
  public:
      /** Full cunstructor. */
    LowRankMatrix(shared_ptr<Matrix<>> A, shared_ptr<Matrix<>> Bt);
    
    /** Matrix-vector-multiply-add: \f$y = y + s A B^\top y \f$. */
    void MultAdd (double s, const BaseVector & x, BaseVector & y) const;

    bool IsComplex() const  { return false; }
    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const;

    virtual int VHeight() const { return m; }
    virtual int VWidth() const { return n; }

    //virtual AutoVector CreateRowVector () const override;
    //virtual AutoVector CreateColVector () const override;
  };

}

#endif

