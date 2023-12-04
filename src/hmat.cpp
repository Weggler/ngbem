#include <solve.hpp>        // everything from ngsolve
#include <cmath>

#include "hmat.hpp"

namespace ngbem
{

  /* This routine computes the geometrical information associated with the
     surface FESpace. It returns two arrays containing 3D-coordinates and cluster
     weights for every boundary dof.
  */
  tuple<Array<double>, Array<Vec<3>>> ComputeClusterData(shared_ptr<FESpace> space,
							 LocalHeap &lh)
  {
    Array<Vec<3>> clcoord;
    Array<double> clweight;
    Array<int> clcounter;
    Array<DofId> mapglob2bnd;
    Array<DofId> mapbnd2glob;

    auto mesh = space->GetMeshAccess();

    // setup global-2-boundary mappings;
    BitArray bnddofs(space->GetNDof());
    bnddofs.Clear();
    for (int i = 0; i < mesh->GetNE(BND); i++)
      {
	Array<DofId> dnums;
	space->GetDofNrs(ElementId(BND, i), dnums);
	for (auto d : dnums)
	  bnddofs.SetBit(d);
      }
    mapglob2bnd.SetSize(space->GetNDof());
    mapglob2bnd = -1;
    int ndof = 0;
    for (int i = 0; i < space->GetNDof(); i++)
      if (bnddofs.Test(i))
	{
	  mapglob2bnd[i] = mapbnd2glob.Size();
	  mapbnd2glob.Append(i);
	  ndof++;
	}

    clcoord.SetSize(ndof);
    clcoord = 0;
    clweight.SetSize(ndof);
    clweight = 0;
    clcounter.SetSize(ndof);
    clcounter = 0;

    // run through all surface elements: clweights and counter
    for (int i = 0; i < mesh->GetNSE(); i++)
      {
	ElementId ei(BND, i);

	Array<DofId> dnumsi;
	space->GetDofNrs(ei, dnumsi); // local dofs
	for (int ii = 0; ii < dnumsi.Size(); ii++)
	  {
	    clweight[mapglob2bnd[dnumsi[ii]]] += mesh->SurfaceElementVolume(i);
	    clcounter[mapglob2bnd[dnumsi[ii]]] += 1;
	  }
      }

    // run through all surface elements: clcoord
    for (int i = 0; i < mesh->GetNSE(); i++)
      {
	HeapReset hr(lh);
	ElementId ei(BND, i);

	/* The distinguished points in an element are its vertices, edge midpoints
	   and element centroid. */
	Array<IntegrationPoint> xhat;
	xhat.Append ( Vec<2> (1, 0) );
	xhat.Append ( Vec<2> (0, 1) );
	xhat.Append ( Vec<2> (0, 0) );
	xhat.Append ( Vec<2> (0.5, 0.5) );
	xhat.Append ( Vec<2> (0, 0.5) );
	xhat.Append ( Vec<2> (0.5, 0) );
	xhat.Append ( Vec<2> (1./3, 1./3) );

	Array<Vec<3>> mipx;
	mipx.SetSize(7);

	// Map the reference points to the physical element
	ElementTransformation &trafoi = mesh->GetTrafo(ei, lh);
	for (int j=0; j<7; j++)
	  {
	    MappedIntegrationPoint<2,3> mip(xhat[j], trafoi);
	    mipx[j] = mip.Point();
	  }

	// run through dofs of the current element
	Array<DofId> dnumsi;
	space->GetDofNrs(ei, dnumsi); // local dofs
	int vertcounter = 0;
	int edgecounter = 0;
	for (int ii = 0; ii < dnumsi.Size(); ii++)
	  {
	    int clc = clcounter[mapglob2bnd[dnumsi[ii]]];
	    // vertex dof is mapped to its vertex
	    if(clc  > 2)
	      {
		clcoord[mapglob2bnd[dnumsi[ii]]] = mipx[vertcounter];
		vertcounter++;
	      }
	    /* Since we do not know how many edge dofs exist, we need to count
	       them first. */
	    if(clc  == 2)
	      {
		edgecounter++;
	      }
	    // bubble dof is mapped to the centroid
	    if(clc  == 1)
	      {
		clcoord[mapglob2bnd[dnumsi[ii]]] = mipx[6];
	      }
	  }

	// In a second loop, map the edge dofs to their edge midpoints
	int edgeoffset = edgecounter/3;
	edgecounter = 0;
	for (int ii = 0; ii < dnumsi.Size(); ii++)
	  {
	    int clc = clcounter[mapglob2bnd[dnumsi[ii]]];
	    if(clc  == 2) // edge dof
	      {
		// The edge dofs of the same edge are mapped to the same point
		for (int offset = 0; offset < edgeoffset; offset++)
		  {
		    clcoord[mapglob2bnd[dnumsi[ii+offset]]] = mipx[3+edgecounter];
		  }

		edgecounter++;
		ii +=edgeoffset-1;
	      }
	  }
      }
    for (int i = 0; i < ndof; i++)
      {
	clweight[i] /= (double) clcounter[i];
	//cout << i << ":" << clweight[i] << "\t" << clcounter[i] << endl;
	//cout << i << ":" << clcoord[i] << endl;
      }


    return tuple {clweight, clcoord};
  }

  /* by S. Rjasanow
     Cluster division algorithm based on principal component analysis.
  */
  void Rja_DivideCluster(Array<double> &G, Array<Vec<3>> &X, int *NClusters, int IClu,
			 Array<struct Cluster> &Clusters, Array<DofId> &Permu, int leafsize)
  {
    // Indices of current cluster
    int NClu = Clusters[IClu].Number;
    int ClBegin = Clusters[IClu].PermuPos;
    int ClEnd = Clusters[IClu].PermuPos + NClu - 1;

    // Moments of cluster
    double Mom[10] = {0.};
    for(int i = ClBegin; i <= ClEnd; i++)
      {
	int j = Permu[i];
	Mom[0] += G[j];
	Mom[1] += G[j] * X[j](0);
	Mom[2] += G[j] * X[j](1);
	Mom[3] += G[j] * X[j](2);
	Mom[4] += G[j] * X[j](0) * X[j](0);
	Mom[5] += G[j] * X[j](0) * X[j](1);
	Mom[6] += G[j] * X[j](0) * X[j](2);
	Mom[7] += G[j] * X[j](1) * X[j](1);
	Mom[8] += G[j] * X[j](1) * X[j](2);
	Mom[9] += G[j] * X[j](2) * X[j](2);
      }

    // Centroid of cluster
    double Centre[3] = {Mom[1] / Mom[0], Mom[2] / Mom[0], Mom[3] / Mom[0]};

    // Covariance matrix of cluster
    double A[9];
    A[0] = Mom[4] - Centre[0] * Mom[1];
    A[1] = Mom[5] - Centre[0] * Mom[2];
    A[2] = Mom[6] - Centre[0] * Mom[3];
    A[3] = A[1];
    A[4] = Mom[7] - Centre[1] * Mom[2];
    A[5] = Mom[8] - Centre[1] * Mom[3];
    A[6] = A[2];
    A[7] = A[5];
    A[8] = Mom[9] - Centre[2] * Mom[3];

    // Get eigenvectors and eigenvalues
    int N = 3, LDA = 3, LWORK = 8, INFO;
    double W[3], WORK[8];
    dsyev_((char *) "V", (char *) "U", &N, A, &LDA, W, WORK, &LWORK, &INFO);

    // Radius and bounding box
    int jMin;
    double XMin[3], XMax[3];
    double Radius = 0., DistMin;
    for(int i = ClBegin; i <= ClEnd; i++)
      {
	int j = Permu[i];

	double Diff[3] = {X[j](0), X[j](1), X[j](2)};
	Diff[0] -= Centre[0];
	Diff[1] -= Centre[1];
	Diff[2] -= Centre[2];

	double Dist = sqrt(fabs(Diff[0] * Diff[0] + Diff[1] * Diff[1] +
				Diff[2] * Diff[2]));

	if(i == ClBegin || Dist < DistMin)
	  {
	    jMin=j;
	    DistMin=Dist;
	  }

	double XT[3];
	XT[0] = A[0] * Diff[0] + A[1] * Diff[1] + A[2] * Diff[2];
	XT[1] = A[3] * Diff[0] + A[4] * Diff[1] + A[5] * Diff[2];
	XT[2] = A[6] * Diff[0] + A[7] * Diff[1] + A[8] * Diff[2];
	Radius = fmax(Radius, sqrt(fabs(XT[0] * XT[0] + XT[1] * XT[1] +
					XT[2] * XT[2])));

	if (i == ClBegin)
	  {
	    XMin[0]=XT[0];
	    XMax[0]=XT[0];
	    XMin[1]=XT[1];
	    XMax[1]=XT[1];
	    XMin[2]=XT[2];
	    XMax[2]=XT[2];
	  }
	else
	  {
	    XMin[0]=fmin(XMin[0], XT[0]);
	    XMax[0]=fmax(XMax[0], XT[0]);
	    XMin[1]=fmin(XMin[1], XT[1]);
	    XMax[1]=fmax(XMax[1], XT[1]);
	    XMin[2]=fmin(XMin[2], XT[2]);
	    XMax[2]=fmax(XMax[2], XT[2]);
	  }
      }

    // Diagonal
    double Diag[3] = {XMax[0] - XMin[0], XMax[1] - XMin[1], XMax[2] - XMin[2]};
    double DiagLength = sqrt(fabs(Diag[0] * Diag[0] + Diag[1] * Diag[1] +
				  Diag[2] * Diag[2]));
    if (DiagLength > 0.0)
      {
	Diag[0] /= DiagLength;
	Diag[1] /= DiagLength;
	Diag[2] /= DiagLength;
      }

    // Copy the data to current cluster
    Clusters[IClu].BasisElement=jMin;
    Clusters[IClu].Radius=Radius;
    Clusters[IClu].DiagLength=DiagLength;
    for (int k = 0; k < 3; k++)
      {
	Clusters[IClu].EVal[k] = W[k];
	Clusters[IClu].XMin[k] = XMin[k];
	Clusters[IClu].XMax[k] = XMax[k];
	Clusters[IClu].Centre[k] = Centre[k];
      }
    for (int k = 0; k < 9; k++)
      Clusters[IClu].EVec[k] = A[k];

    // If the size cluster is smaller than leafsize, it becomes a leaf (no children)
    if (NClu <= leafsize)
      {
	Clusters[IClu].Child1 =- 1;
	Clusters[IClu].Child2 =- 1;
	return;
      }
    // Otherwise we divide the cluster into two children
    else
      {
	int Num1 = 0, Num2 = 0;
	int ICheck = ClBegin;
	int IEnd = ClEnd;
	double Crit = Centre[0] * A[6] + Centre[1] * A[7] + Centre[2] * A[8];

	/* Split along the eigenvector of the largest eigenvalue of the
	   covariance matrix. */
	for (int i = 0; i < NClu; i++)
	  {
	    int j = Permu[ICheck];
	    if (X[j](0) * A[6] + X[j](1) * A[7] + X[j](2) * A[8] >= Crit)
	      {
		Num1 += 1;
		ICheck += 1;
	      }
	    else
	      {
		Num2 += 1;
		Permu[ICheck] = Permu[IEnd];
		Permu[IEnd] = j;
		IEnd -= 1;
	      }
	  }

	// Non-separable cluster
	if (Num1 == 0 || Num2 == 0)
	  {
	    Clusters[IClu].Child1 = -1;
	    Clusters[IClu].Child2 = -1;
	    return;
	  }
	// Separable cluster with two children
	else
	  {
	    // Provide hierarchical information
	    struct Cluster cl1, cl2;
	    int IClu1 = *NClusters;
	    Clusters[IClu].Child1 = IClu1;
	    cl1.Level = Clusters[IClu].Level + 1;
	    cl1.Parent = IClu;
	    cl1.Number = Num1;
	    cl1.PermuPos = ClBegin;
	    Clusters.Append(cl1);

	    int IClu2 = *NClusters + 1;
	    Clusters[IClu].Child2 = IClu2;
	    cl2.Level = Clusters[IClu].Level + 1;
	    cl2.Parent = IClu;
	    cl2.Number = Num2;
	    cl2.PermuPos = ClBegin + Num1;
	    Clusters.Append(cl2);

	    *NClusters += 2;
	  }
      }
  }

  /* by S. Rjasanow
     Generate a cluster tree from weights and points. Returns the number of clusters.
  */
  int Rja_ClusterTree(int n_dof, Array<double> &G, Array<Vec<3>> &X,
		      Array<struct Cluster> &Clusters, Array<DofId> &Permu, int leafsize)
  {
    /* Initial cluster */
    int NClusters = 1;
    int IClu = 0;

    /* Initial permutation */
    for (int i = 0; i < n_dof; i++)
      Permu[i]=i;

    /* Initial data */
    struct Cluster cl1;
    cl1.Level = 0;
    cl1.Parent = -1;
    cl1.Number = n_dof;
    cl1.PermuPos = 0;
    Clusters.Append(cl1);
    
    /* Devide cluster */
    do
      {
	Rja_DivideCluster(G, X, &NClusters, IClu, Clusters, Permu, leafsize);
	IClu += 1;
      }
    while (IClu < NClusters);

    return NClusters;
  }

  ClusterTree :: ClusterTree(shared_ptr<FESpace> _space, int leafsize)
    : space(_space)
  {
    // Get points and weights for clustering
    LocalHeap locheap(1000000);
    auto [G, X] = ComputeClusterData(space, locheap);

    // Number of dofs, i.e. length of index set
    int n_dof = G.Size();

    // Mapping from FE ordering of dofs to cluster ordering of dofs
    mapbnd2cluster.SetSize(n_dof);

    // Generate a cluster tree
    n_cluster = Rja_ClusterTree(n_dof, G, X, arr_clusters, mapbnd2cluster, leafsize);

    cout << "number of clusters: " << n_cluster << endl;
    cout << "bnd2cluster: " << mapbnd2cluster << endl;
  }

}
