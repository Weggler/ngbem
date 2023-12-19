/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                16.03.2009                                   *
\*****************************************************************************/
/*
  Last change 16.03.2009
*/

/*
  Includes for all functions 
*/

# include "Includes.h"

long AdmissibleBlocks(Cluster ClCand,double *Points,double *Rows, long *aNRows)
{
/*
  Local variables
*/
  double Dist[3],DistNorm,norm,R;
  long NRows;

  double alpha,beta,gamma;

  alpha=0.1;
  beta=0.8;
  gamma=1.0;

  NRows=*aNRows;
    
/*
  Distance between the clusters
*/
  R=ClCand.Radius;

  IEnd=NRows;
  NCheck=0;

  do
    {
      cblas_dcopy(3,ClCand.Centre,1,Dist,1);
      cblas_daxpy(3,-1.0,Points[3*Rows[NCheck]],1,Dist,1);
      DistNorm=cblas_dnrm2(3,Dist,1)-R;
/*
  1. Criterium
*/
      if (DistNorm >= alpha*Kappa*R*R+beta*R)
	NCheck+=1;
      else
	{
	  ih=Rows[IEnd];
	  Rows[IEnd]=Rows[NCheck];
	  Rows[NCheck]=ih;
	  IEnd-=1;
	}
    }
  while (IEnd > NCheck);

  if (NCheck == 0)
    return 0;
/*
  Some candidates found
*/
  else
    {
      NAdPoints=NCheck;
/*
  Projections
*/
      AdPoints=(double *)malloc(3*NAdPoints*sizeof(double));
      assert(AdPoints != NULL);

      for (i=0; i < NAdPoints; i++)
	{
	  cblas_dcopy(3,Points+3*Rows[i],1,AdPoints+3*i,1);
	  norm=1.0/cblas_dnrm2(3,AdPoints+3*i,1);
	  cblas_dscal(3,norm,AdPoints+3*i,1);
	}
/*
  Clustering
*/
      GClusterTree();
    }      
}
