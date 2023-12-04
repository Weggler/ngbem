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

    /** The clustering produces a different ordering of boundary dofs. */
    Array<DofId> mapbnd2cluster;

    /** Number of clusters in the tree. */
    int n_cluster;

    /** Array of length #n_cluster containing the clusters. */
    Array<struct Cluster> arr_clusters;

  public:
    ClusterTree(shared_ptr<FESpace> space, int leafsize);
  };

  // class HMatrix
  // {
  //   shared_ptr<ClusterTree> row_cluster;
  //   shared_ptr<ClusterTree> col_cluster;
    
  // public:
  //   HMatrix(shared_ptr<ClusterTree> row_cluster, shared_ptr<ClusterTree> col_cluster, double eta);
  // };

}



#endif

