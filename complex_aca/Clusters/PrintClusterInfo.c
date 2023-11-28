/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                02.01.2009                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void PrintClusterInfo(long NClusters,Cluster *Clusters,long *Permu,long NPairs,Pair *Pairs)
{
/*
  Local variables
*/
    FILE *ResultFile;
    long i;

    ResultFile=fopen("surface.cls","w");
/*
  Number of elements
*/
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"# number of elements\n");
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"%ld\n",NumElements);
/*
  Permutation
*/
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"# permutation\n");
    fprintf(ResultFile,"#\n");
    for (i=0; i < NumElements; i++)
      fprintf(ResultFile,"%ld %ld\n",i+1,Permu[i]+1);
/*
  Clusters
*/
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"# number of clusters\n");
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"%ld\n",NClusters);
    for (i=0; i < NClusters; i++)
	fprintf(ResultFile,"%ld %ld %ld %ld\n",i+1,Clusters[i].Level,Clusters[i].Number,Clusters[i].PermuPos+1);
/*
  Pairs
*/
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"# number of pairs\n");
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"%ld\n",NPairs);
    for (i=0; i < NPairs; i++)
      fprintf(ResultFile,"%ld %ld %ld %ld\n",i+1,Pairs[i].Clu1+1,Pairs[i].Clu2+1,(long)Pairs[i].Type);

    fclose(ResultFile);
}
