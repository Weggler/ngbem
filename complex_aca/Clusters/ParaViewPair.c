/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                24.03.2010                                   *
\*****************************************************************************/
/*
  Last change 24.03.2010
*/

/*
  Includes for all functions 
*/

# include "Includes.h"

void ParaViewPair(long PairNumber,Pair *Pairs,Cluster *Clusters,long *Permu)
{
/*
  Local variables
*/
  char ResultName[4096];
  FILE *ResultFile=0;
  long i,NClu,Pos;
  long *Type;
/*
  Memory allocation
*/
  Type=(long *)malloc(NFeko*sizeof(long));
  assert(Type != NULL);
  
  sprintf(ResultName,"ParaVeiwCluster%ld.dat",PairNumber);
  ResultFile=fopen(ResultName,"w");

  fprintf(ResultFile,"BEG_ELM_VARS\n");
  fprintf(ResultFile,"Type\n");

  lset(NFeko,1,Type,1);
/*
  First cluster
*/
  NClu=Clusters[Pairs[PairNumber].Clu1].Number; 
  Pos=Clusters[Pairs[PairNumber].Clu1].PermuPos;

  printf("\n 1st: %ld %ld\n",NClu,Pos);

  for (i=0; i < NClu; i++)
    Type[Permu[Pos+i]]=2;
/*
  Second cluster
*/
  NClu=Clusters[Pairs[PairNumber].Clu2].Number; 
  Pos=Clusters[Pairs[PairNumber].Clu2].PermuPos;

  printf(" 2nd: %ld %ld\n",NClu,Pos);

  for (i=0; i < NClu; i++)
    Type[Permu[Pos+i]]=2;

  for (i=0; i < NFeko; i++)
    fprintf(ResultFile,"%ld\n",Type[i]);

  fclose(ResultFile);

  free(Type);
}
