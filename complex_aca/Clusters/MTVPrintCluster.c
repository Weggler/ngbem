/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                02.04.2009                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void MTVPrintCluster(long N,double *X,Cluster *Clusters,long *Permu,Pair ClPair)
{
/*
  Local variables
*/
    FILE *ResultFile;
    long i,j,NClu;

    ResultFile=fopen("ClusterPair.mtv","w");
/*
  Header
*/
    fprintf(ResultFile,"$ data = curve3d\n");
    fprintf(ResultFile,"%% mt=1 mc=1 lt=0 equalscale\n");
/*
  All points
*/
    for (i=0; i < N; i++)
      fprintf(ResultFile,"%14.8e %14.8e %14.8e\n",X[3*i],X[3*i+1],X[3*i+2]);
/*
  Second header
*/
    fprintf(ResultFile,"\n");
    fprintf(ResultFile,"%% mt=1 mc=4 lt=0 equalscale\n");
/*
  First cluster
*/
    NClu=Clusters[ClPair.Clu1].Number;
    for (i=0; i < NClu; i++)
      {
	j=Permu[Clusters[ClPair.Clu1].PermuPos+i];
        fprintf(ResultFile,"%14.8e %14.8e %14.8e\n",X[3*j],X[3*j+1],X[3*j+2]);
      }
/*
  Second cluster
*/
    NClu=Clusters[ClPair.Clu2].Number;
    for (i=0; i < NClu; i++)
      {
	j=Permu[Clusters[ClPair.Clu2].PermuPos+i];
        fprintf(ResultFile,"%14.8e %14.8e %14.8e\n",X[3*j],X[3*j+1],X[3*j+2]);
      }

    fclose(ResultFile);
}
