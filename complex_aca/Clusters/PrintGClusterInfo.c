/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                06.04.2009                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void PrintGClusterInfo(long N,long NClusters,Cluster *Clusters,long *Permu)
{
/*
  Local variables
*/
    FILE *ResultFile;
    long i;

    ResultFile=fopen("GSurface.cls","w");
/*
  Number of elements
*/
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"# number of elements\n");
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"%ld\n",N);
/*
  Permutation
*/
    fprintf(ResultFile,"#\n");
    fprintf(ResultFile,"# permutation\n");
    fprintf(ResultFile,"#\n");
    for (i=0; i < N; i++)
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

    fclose(ResultFile);
}
