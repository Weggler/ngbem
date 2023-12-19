/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                12.11.2006                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void GetHMatrix(long N,double *X,double *G,HMatrix *HMat,
                double complex (*ElementName)(long,long),
                ACAParameters *aParamACA)
{
/*
  Local variables
*/
  long NClusters,NPairs;
  Cluster *Clusters;
  Pair *Pairs;
  long *Permu,NMax,MMax;
  ACAParameters ParamACA;
/* 
  Parameters
*/
  ParamACA=*aParamACA;
/* 
  Construct cluster tree for the points
*/
  if (ParamACA.Info == 1)
    {
      printf("Construct cluster tree ...\n"); 
      fflush(stdout);
    }

  NClusters=ClusterTree(N,G,X,&Clusters,&Permu,ParamACA);

  if (ParamACA.Info == 1)
    {
      printf("  Done with %ld clusters\n\n",NClusters);
      fflush(stdout);

      printf("Construct cluster pairs ...\n"); 
      fflush(stdout);
    }

  NPairs=SymmClusterPairs(Clusters,&Pairs,&NMax,&MMax,ParamACA);

  ParamACA.NMax=NMax;
  ParamACA.MMax=MMax;

  if (ParamACA.Info == 1)
    {
      printf("  Done with %ld cluster pairs\n",NPairs);
      printf("  Maximal dimension NMax = %ld, MMax = %ld\n\n",NMax,MMax);
      fflush(stdout);
 
      printf("Generate H-Matrix ...\n"); 
      fflush(stdout);
    }

  GenHMatrix(N,N,0,0,
             Permu,Permu,
	     Clusters,Clusters,
	     NPairs,Pairs,HMat,
	     ElementName,&ParamACA);
  
  *aParamACA=ParamACA;
/*
  Free memory
*/
  free(Clusters); free(Pairs);
}
