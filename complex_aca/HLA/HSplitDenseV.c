/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                10.12.2010                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HSplitDenseV(CBlock *A,long M11,CBlock *A11,CBlock *A12)
{
  HCheckM(A);
  HCheckM(A11);
  HCheckM(A12);
/*
  Local variables
*/
  long NA,MA,N11,N12,M12;

  assert((*A).Type == Dense);

  NA=(*A).NBlock;
  MA=(*A).MBlock;
/*
  Block 1,1
*/
  N11=NA;

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

  cblas_zcopy(N11*M11,(*A).Data,1,(*A11).Data,1);

  (*A11).Rank=0;
  (*A11).UData=NULL;
  (*A11).VData=NULL;
/*
  Block 1,2
*/
  N12=NA;
  M12=MA-M11;
  assert(M12 >= 1);

  (*A12).IBlock=(*A).IBlock;
  (*A12).JBlock=(*A).JBlock+M11;
  (*A12).NBlock=N12;
  (*A12).MBlock=M12;

  (*A12).Type=Dense;

  (*A12).A11=NULL;
  (*A12).A21=NULL;
  (*A12).A12=NULL;
  (*A12).A22=NULL;

  (*A12).Data=(double complex *)malloc(N12*M12*sizeof(double complex));
  assert((*A12).Data != NULL);

  cblas_zcopy(N12*M12,(*A).Data+N11*M11,1,(*A12).Data,1);

  (*A12).Rank=0;
  (*A12).UData=NULL;
  (*A12).VData=NULL;

  HCheckM(A);
  HCheckM(A11);
  HCheckM(A12);
}
