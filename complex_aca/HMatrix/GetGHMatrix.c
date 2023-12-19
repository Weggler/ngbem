/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                06.04.2009                                   *
\*****************************************************************************/
/*
  Includes for the main function 
*/

# include "../CInclude/Includes.h"

void GetGHMatrix(long N,double *X,double *G,GHMatrix *GHMat,
                 double complex (*ElementName)(long,long),
                 ACAParameters *aParamACA)
{
/*
  Local variables
*/
  long NClusters,NBlocks;
  Cluster *Clusters;
  GCBlock *Blocks;
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
  PrintGClusterInfo(N,NClusters,Clusters,Permu);

  if (ParamACA.Info == 1)
    {
      printf("  Done with %ld clusters\n\n",NClusters);
      fflush(stdout);

      printf("Construct blocks ...\n"); 
      fflush(stdout);
    }

  NBlocks=ClusterPoint(N,Clusters,Permu,X,
                       &Blocks,&NMax,&MMax,ParamACA);

  (*GHMat).NMax=NMax;
  (*GHMat).MMax=MMax;

  if (ParamACA.Info == 1)
    {
      printf("  Done with %ld blocks\n",NBlocks);
      printf("  Maximal dimension NMax = %ld, MMax = %ld\n\n",NMax,MMax);
      fflush(stdout);

 
      printf("Generate H-Matrix ...\n"); 
      fflush(stdout);
    }

  GenGHMatrix(N,N,
	      NBlocks,Blocks,GHMat,
	      ElementName,&ParamACA);

  *aParamACA=ParamACA;
}
