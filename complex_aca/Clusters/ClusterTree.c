//                                                                  
// S. Rjasanow: Adaptive Cross Approximation                 
//
//------------------------------------------------------------------
//                                                                  
//  routine name  -  ClusterTree
//                                                                  
//------------------------------------------------------------------
//                                                                  
//  last revision -  Aug 12
//
//   purpose      - generate a cluster from geometrical data
//
//   in:
//              N - number of points to be clustered
//              G - cluster weights 
//              X - cluster points 
//       ParamACA - ACA parameters used to generate cluster
//
//   out:
//      aClusters - cluster tree 
//         aPermu - permutation in X,G ... (?)
//                                                                  
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

long ClusterTree(long N,double *G,double *X,
                 Cluster **aClusters, long **aPermu,
                 ACAParameters ParamACA)
{
  /* Local variables */
  long i,IClu,NClusters,*Permu;
  Cluster *Clusters;

  assert(N >= 1);

  Permu=(long *)malloc(N*sizeof(long));
  assert(Permu != NULL);

  /* Initial cluster */
  NClusters=1;
  IClu=0;

  Clusters=(Cluster *)malloc(NClusters*sizeof(Cluster));
  assert(Clusters != NULL);

  /* Initial permutation */
  for (i=0; i < N; i++)
    Permu[i]=i;

  /* Initial data */
  Clusters[IClu].Level = 0;
  Clusters[IClu].Father = -1;
  Clusters[IClu].Number = N;
  Clusters[IClu].PermuPos = 0;

  /* Devide cluster */
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
