/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                11.03.2011                                   *
\*****************************************************************************/
/*
  Last change 05.04.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

double HNormM(CBlock *A)
{
  double complex FroNorm,CU,CV;
  long k,l;

  if((*A).Type == Null)
    return 0.0;
  else if((*A).Type  == Hierarchical)
    return sqrt(pow(HNormM((*A).A11),2)+pow(HNormM((*A).A21),2)+pow(HNormM((*A).A12),2)+pow(HNormM((*A).A22),2));
  else if ((*A).Type == Dense)
    {
      cblas_zdotc_sub((*A).NBlock*(*A).MBlock,(*A).Data,1,(*A).Data,1,&FroNorm);
      return sqrt(creal(FroNorm));
    }
  else if ((*A).Type == Admissible)
    {
      FroNorm=0.0+I*0.0;

      for (l=0; l < (*A).Rank; l++)
	for (k=0; k < (*A).Rank; k++)
	  {
	    cblas_zdotc_sub((*A).NBlock,(*A).UData+k*(*A).NBlock,1,(*A).UData+l*(*A).NBlock,1,&CU);
	    cblas_zdotc_sub((*A).MBlock,(*A).VData+l*(*A).MBlock,1,(*A).VData+k*(*A).MBlock,1,&CV);
	    FroNorm+=CU*CV;
	  }
      return sqrt(creal(FroNorm));
    }
  return 0.0;
}
