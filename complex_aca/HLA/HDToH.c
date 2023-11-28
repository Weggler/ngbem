/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                29.11.2010                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/
/*
  Includes for all functions 
*/

# include "Includes.h"

void HDToH(CBlock *A,long N11,long M11)
{
  HCheckM(A);
/*
  Local variabls
*/
  long NA,MA,N21,M21,N12,M12,N22,M22;
  long j;

  assert((*A).Type == Dense);

  NA=(*A).NBlock;
  MA=(*A).MBlock;

  N21=NA-N11;
  M21=M11;

  N12=N11;
  M12=MA-M11;

  N22=NA-N11;
  M22=MA-M11;

  (*A).Type=Hierarchical;
/*
  Block 1,1
*/
  (*A).A11=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A11 != NULL);
  HIniM((*A).A11);

  (*(*A).A11).IBlock=(*A).IBlock;
  (*(*A).A11).JBlock=(*A).JBlock;
  (*(*A).A11).NBlock=N11;
  (*(*A).A11).MBlock=M11;

  (*(*A).A11).Type=Dense;

  (*(*A).A11).A11=NULL;
  (*(*A).A11).A21=NULL;
  (*(*A).A11).A12=NULL;
  (*(*A).A11).A22=NULL;

  (*(*A).A11).Data=(double complex *)malloc(N11*M11*sizeof(double complex));
  assert((*(*A).A11).Data != NULL);

  for (j=0; j < M11; j++)
    cblas_zcopy(N11,(*A).Data+j*NA,1,(*(*A).A11).Data+j*N11,1);

  (*(*A).A11).Rank=0;
  (*(*A).A11).UData=NULL;
  (*(*A).A11).VData=NULL;
/*
  Block 2,1
*/
  (*A).A21=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A21 != NULL);
  HIniM((*A).A21);

  (*(*A).A21).IBlock=(*A).IBlock+N11;
  (*(*A).A21).JBlock=(*A).JBlock;
  (*(*A).A21).NBlock=N21;
  (*(*A).A21).MBlock=M21;

  (*(*A).A21).Type=Dense;

  (*(*A).A21).A11=NULL;
  (*(*A).A21).A21=NULL;
  (*(*A).A21).A12=NULL;
  (*(*A).A21).A22=NULL;

  (*(*A).A21).Data=(double complex *)malloc(N21*M21*sizeof(double complex));
  assert((*(*A).A21).Data != NULL);

  for (j=0; j < M21; j++)
    cblas_zcopy(N21,(*A).Data+j*NA+N11,1,(*(*A).A11).Data+j*N21,1);

  (*(*A).A21).Rank=0;
  (*(*A).A21).UData=NULL;
  (*(*A).A21).VData=NULL;
/*
  Block 1,2
*/
  (*A).A12=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A12 != NULL);
  HIniM((*A).A12);

  (*(*A).A12).IBlock=(*A).IBlock;
  (*(*A).A12).JBlock=(*A).JBlock+M11;
  (*(*A).A12).NBlock=N12;
  (*(*A).A12).MBlock=M12;

  (*(*A).A12).Type=Dense;

  (*(*A).A12).A11=NULL;
  (*(*A).A12).A21=NULL;
  (*(*A).A12).A12=NULL;
  (*(*A).A12).A22=NULL;

  (*(*A).A12).Data=(double complex *)malloc(N12*M12*sizeof(double complex));
  assert((*(*A).A12).Data != NULL);

  for (j=0; j < M12; j++)
    cblas_zcopy(N12,(*A).Data+(j+M11)*NA,1,(*(*A).A11).Data+j*N12,1);

  (*(*A).A12).Rank=0;
  (*(*A).A12).UData=NULL;
  (*(*A).A12).VData=NULL;
/*
  Block 2,2
*/
  (*A).A22=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A22 != NULL);
  HIniM((*A).A22);

  (*(*A).A22).IBlock=(*A).IBlock+N11;
  (*(*A).A22).JBlock=(*A).JBlock+M11;
  (*(*A).A22).NBlock=N22;
  (*(*A).A22).MBlock=M22;

  (*(*A).A22).Type=Dense;

  (*(*A).A22).A11=NULL;
  (*(*A).A22).A21=NULL;
  (*(*A).A22).A12=NULL;
  (*(*A).A22).A22=NULL;

  (*(*A).A22).Data=(double complex *)malloc(N22*M22*sizeof(double complex));
  assert((*(*A).A22).Data != NULL);

  for (j=0; j < M22; j++)
    cblas_zcopy(N22,(*A).Data+(j+M11)*NA+N11,1,(*(*A).A11).Data+j*N22,1);

  (*(*A).A22).Rank=0;
  (*(*A).A22).UData=NULL;
  (*(*A).A22).VData=NULL;
 
  free((*A).Data);

  (*A).Data=NULL;

  (*A).Rank=0;
  (*A).UData=NULL;
  (*A).VData=NULL;

  HCheckM(A);
}
