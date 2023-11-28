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

long AdmissiblePair(long IClu1,long IClu2,Cluster *Clusters1,Cluster *Clusters2,ACAParameters ParamACA)
{
/*
  Local variables
*/
    double Dist[3],DistNorm,R1,R2;
/*
  Distance between the clusters
*/
    R1=Clusters1[IClu1].Radius;
    R2=Clusters2[IClu2].Radius;

    cblas_dcopy(3,Clusters1[IClu1].Centre,1,Dist,1);
    cblas_daxpy(3,-1.0,Clusters2[IClu2].Centre,1,Dist,1);
    DistNorm=cblas_dnrm2(3,Dist,1);
    DistNorm=max(0.0,DistNorm-R1-R2);
//    printf("%25.16lf %25.16lf %25.16lf\n", R1,R2, DistNorm);
/*
  Criterium
*/
    if (min(R1,R2) <= ParamACA.ACAEta*DistNorm) 
/*    if (R1+R2 <= ParamACA.ACAEta*DistNorm) */
      return 1;
    else
      return 0;
}
