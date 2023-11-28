/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                27.12.2008                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void PrintCluster(long IClu,Cluster *Clusters)
{
    printf("\n");
    printf("Begin cluster print :\n");
    printf("  Cluster number     = %ld\n",IClu);
    printf("  Cluster level      = %ld\n",Clusters[IClu].Level);
    printf("  Cluster father     = %ld\n",Clusters[IClu].Father);
    printf("  Cluster son 1      = %ld\n",Clusters[IClu].Son1);
    printf("  Cluster son 2      = %ld\n",Clusters[IClu].Son2);
    printf("  Number of elements = %ld\n",Clusters[IClu].Number);
    printf("  Position in Permu  = %ld\n",Clusters[IClu].PermuPos);
    printf("  BBMin    = %16.8e %16.8e %16.8e\n",Clusters[IClu].XMin[0],Clusters[IClu].XMin[1],Clusters[IClu].XMin[2]);
    printf("  BBMax    = %16.8e %16.8e %16.8e\n",Clusters[IClu].XMax[0],Clusters[IClu].XMax[1],Clusters[IClu].XMax[2]);
    printf("  Centre   = %16.8e %16.8e %16.8e\n",Clusters[IClu].Centre[0],Clusters[IClu].Centre[1],Clusters[IClu].Centre[2]);
    printf("  Diagonal = %16.8e\n",Clusters[IClu].DiagLength);
    printf("End cluster print\n");
    printf("\n");
}
