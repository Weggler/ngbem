/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                27.03.2009                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

long KappaAdmissiblePair(long IClu1,long IClu2,
                         Cluster *Clusters1,Cluster *Clusters2,
                         ACAParameters ParamACA)
{
/*
  Local variables
*/
  double BB[24],C[3],XMin[3],XMax[3],EV[9];
  double Dist[3],Dist1[3],DistNorm,R1,R2,Phi;
  long i;
/*
  Distance between the clusters
*/
  R1=Clusters1[IClu1].Radius;
  R2=Clusters2[IClu2].Radius;

  cblas_dcopy(3,Clusters1[IClu1].Centre,1,Dist,1);
  cblas_daxpy(3,-1.0,Clusters2[IClu2].Centre,1,Dist,1);
  DistNorm=cblas_dnrm2(3,Dist,1);
  DistNorm=max(0.0,DistNorm-R1-R2);

  if (DistNorm <= ParamACA.ACAAlpha*R1+ParamACA.ACABeta*Kappa*R1*R1)
    return 0;
  else
    {
/*
  Bounding box, 8 corners
*/
      cblas_dcopy(3,Clusters1[IClu2].Centre,1,C,1);
      cblas_dcopy(3,Clusters1[IClu2].XMin,1,XMin,1);
      cblas_dcopy(3,Clusters1[IClu2].XMax,1,XMax,1);
      cblas_dcopy(9,Clusters1[IClu2].EVec,1,EV,1);

      cblas_dcopy(3,C,1,BB+0,1);
      cblas_daxpy(3,XMin[0],EV+0,1,BB+0,1);
      cblas_daxpy(3,XMin[1],EV+3,1,BB+0,1);
      cblas_daxpy(3,XMin[2],EV+6,1,BB+0,1);

      cblas_dcopy(3,C,1,BB+3,1);
      cblas_daxpy(3,XMin[0],EV+0,1,BB+3,1);
      cblas_daxpy(3,XMax[1],EV+3,1,BB+3,1);
      cblas_daxpy(3,XMin[2],EV+6,1,BB+3,1);

      cblas_dcopy(3,C,1,BB+6,1);
      cblas_daxpy(3,XMax[0],EV+0,1,BB+6,1);
      cblas_daxpy(3,XMin[1],EV+3,1,BB+6,1);
      cblas_daxpy(3,XMin[2],EV+6,1,BB+6,1);

      cblas_dcopy(3,C,1,BB+9,1);
      cblas_daxpy(3,XMax[0],EV+0,1,BB+9,1);
      cblas_daxpy(3,XMax[1],EV+3,1,BB+9,1);
      cblas_daxpy(3,XMin[2],EV+6,1,BB+9,1);

      cblas_dcopy(3,C,1,BB+12,1);
      cblas_daxpy(3,XMin[0],EV+0,1,BB+12,1);
      cblas_daxpy(3,XMin[1],EV+3,1,BB+12,1);
      cblas_daxpy(3,XMax[2],EV+6,1,BB+12,1);

      cblas_dcopy(3,C,1,BB+15,1);
      cblas_daxpy(3,XMin[0],EV+0,1,BB+15,1);
      cblas_daxpy(3,XMax[1],EV+3,1,BB+15,1);
      cblas_daxpy(3,XMax[2],EV+6,1,BB+15,1);

      cblas_dcopy(3,C,1,BB+18,1);
      cblas_daxpy(3,XMax[0],EV+0,1,BB+18,1);
      cblas_daxpy(3,XMin[1],EV+3,1,BB+18,1);
      cblas_daxpy(3,XMax[2],EV+6,1,BB+18,1);

      cblas_dcopy(3,C,1,BB+21,1);
      cblas_daxpy(3,XMax[0],EV+0,1,BB+21,1);
      cblas_daxpy(3,XMax[1],EV+3,1,BB+21,1);
      cblas_daxpy(3,XMax[2],EV+6,1,BB+21,1);
/*
  Maximal angle
*/
      Phi=0.0;

      for (i=0; i < 8; i++)
	{
	  cblas_dcopy(3,Clusters1[IClu1].Centre,1,Dist1,1);
	  cblas_daxpy(3,-1.0,BB+3*i,1,Dist1,1);

	  Phi=max(Phi,acos(cblas_ddot(3,Dist,1,Dist1,1)/
			   (cblas_dnrm2(3,Dist,1)*cblas_dnrm2(3,Dist1,1))));
	}
/*
  Final decision
*/
      if (Phi <= M_PI/(1.0+ParamACA.ACAGamma*R1*Kappa))
	return 1;
      else
	return 0;
    }
}
