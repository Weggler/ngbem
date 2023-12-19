/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                20.03.2009                                   *
\*****************************************************************************/
/*
  Last change 16.03.2009
*/

/*
  Includes for all functions 
*/

# include "Includes.h"

long GClusterTree(long N,double *X,Cluster **aClusters, long **aPermu,ACAParameters ParamACA)
{
/*
  Local variables
*/
  long i,IClu,NClusters,*Permu;
  Cluster *Clusters;

  if (N < 1)
    {
      printf("Wrong number of points in ClusterTree, N=%ld\n",N);
      exit(0);
    }

  Permu=(long *)malloc(N*sizeof(long));
  assert(Permu != NULL);
/*
  Initial cluster
*/
  NClusters=1;
  IClu=0;

  Clusters=(Cluster *)malloc(NClusters*sizeof(Cluster));
  assert(Clusters != NULL);
/*
  Initial permutation
*/
  for (i=0; i < N; i++)
    Permu[i]=i;
/*
  Initial data
*/
  Clusters[IClu].Level=0;
  Clusters[IClu].Father=-1;
  Clusters[IClu].Number=N;
  Clusters[IClu].PermuPos=0;
/*
  Devide cluster
*/
  do 
    {
      DevideCluster(G,X,&NClusters,IClu,&Clusters,Permu,ParamACA);
      IClu+=1;
    }
  while (IClu < NClusters);

  *aClusters=Clusters;
  *aPermu=Permu;

  return NClusters;   
}
