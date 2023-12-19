/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                10.12.2010                                   *
\*****************************************************************************/
/*
  Last change 10.12.2010
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HSplitDenseH(CBlock *A,long N11,CBlock *A11,CBlock *A21)
{
  HCheckM(A);
  HCheckM(A11);
  HCheckM(A21);
/*
  Local variables
*/
  long NA,MA,M11,N21,M21;
  long j;

  assert((*A).Type == Dense);

  NA=(*A).NBlock;
  MA=(*A).MBlock;
/*
  Block 1,1
*/
  M11=MA;

  (*A11).IBlock=(*A).IBlock;
  (*A11).JBlock=(*A).JBlock;
  (*A11).NBlock=N11;
  (*A11).MBlock=M11;

  (*A11).Type=Dense;

  (*A11).A11=NULL;
  (*A11).A21=NULL;
  (*A11).A12=NULL;
  (*A11).A22=NULL;

  (*A11).Data=(double complex *)malloc(N11*M11*sizeof(double complex));
  assert((*A11).Data != NULL);

  for (j=0; j < M11; j++)
    cblas_zcopy(N11,(*A).Data+j*NA,1,(*A11).Data+j*N11,1);

  (*A11).Rank=0;
  (*A11).UData=NULL;
  (*A11).VData=NULL;
/*
  Block 2,1
*/
  N21=NA-N11;
  M21=MA;
  assert(N21 >= 1);

  (*A21).IBlock=(*A).IBlock;
  (*A21).JBlock=(*A).JBlock+M11;
  (*A21).NBlock=N21;
  (*A21).MBlock=M21;

  (*A21).Type=Dense;

  (*A21).A11=NULL;
  (*A21).A21=NULL;
  (*A21).A12=NULL;
  (*A21).A22=NULL;

  (*A21).Data=(double complex *)malloc(N21*M21*sizeof(double complex));
  assert((*A21).Data != NULL);

  for (j=0; j < M21; j++)
    cblas_zcopy(N21,(*A).Data+j*NA+N11,1,(*A21).Data+j*N21,1);

  (*A21).Rank=0;
  (*A21).UData=NULL;
  (*A21).VData=NULL;

  HCheckM(A);
  HCheckM(A11);
  HCheckM(A21);
}
