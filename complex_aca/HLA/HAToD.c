/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                08.12.2010                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HAToD(CBlock *A)
{
  HCheckM(A);
/*
  Local variables
*/
  long NA,MA;
  double complex C0,C1;

  assert((*A).Type == Admissible);
/*
  Initialisation
*/
  C0=0.0+0.0*I;
  C1=1.0+0.0*I;

  (*A).Type=Dense;

  NA=(*A).NBlock;
  MA=(*A).MBlock;

  (*A).Data=(double complex *)malloc(NA*MA*sizeof(double complex));
  assert((*A).Data != NULL);

  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,NA,MA,(*A).Rank,&C1,(*A).UData,NA,(*A).VData,MA,&C0,(*A).Data,NA);

  free((*A).UData); free((*A).VData);

  (*A).Rank=0;
  (*A).UData=NULL;
  (*A).VData=NULL;

  HCheckM(A);
}
