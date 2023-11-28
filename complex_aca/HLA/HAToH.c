/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                01.12.2010                                   *
\*****************************************************************************/
/*
  Last change 06.12.2010
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HAToH(CBlock *A,long N11,long M11)
{
  HCheckM(A);
/*
  Local variables
*/
  long NA,MA,N21,M21,N12,M12,N22,M22,Rank;
  long j;

  assert((*A).Type == Admissible);

  NA=(*A).NBlock;
  MA=(*A).MBlock;

  Rank=(*A).Rank;

  N21=NA-N11;
  assert(N21 >= 1);
  M21=M11;

  N12=N11;
  M12=MA-M11;
  assert(M12 >= 1);

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

  (*(*A).A11).Type=Admissible;

  (*(*A).A11).A11=NULL;
  (*(*A).A11).A21=NULL;
  (*(*A).A11).A12=NULL;
  (*(*A).A11).A22=NULL;

  (*(*A).A11).Data=NULL;

  (*(*A).A11).Rank=Rank;

  (*(*A).A11).UData=(double complex *)malloc(N11*Rank*sizeof(double complex));
  assert((*(*A).A11).UData != NULL);
  (*(*A).A11).VData=(double complex *)malloc(M11*Rank*sizeof(double complex));
  assert((*(*A).A11).UData != NULL);

  for (j=0; j < Rank; j++)
    {
      cblas_zcopy(N11,(*A).UData+j*NA,1,(*(*A).A11).UData+j*N11,1);
      cblas_zcopy(M11,(*A).VData+j*MA,1,(*(*A).A11).VData+j*M11,1);
    }
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

  (*(*A).A21).Type=Admissible;

  (*(*A).A21).A11=NULL;
  (*(*A).A21).A21=NULL;
  (*(*A).A21).A12=NULL;
  (*(*A).A21).A22=NULL;

  (*(*A).A21).Data=NULL;

  (*(*A).A21).Rank=Rank;

  (*(*A).A21).UData=(double complex *)malloc(N21*Rank*sizeof(double complex));
  assert((*(*A).A21).UData != NULL);
  (*(*A).A21).VData=(double complex *)malloc(M21*Rank*sizeof(double complex));
  assert((*(*A).A21).UData != NULL);

  for (j=0; j < Rank; j++)
    {
      cblas_zcopy(N21,(*A).UData+j*NA+N11,1,(*(*A).A21).UData+j*N21,1);
      cblas_zcopy(M21,(*A).VData+j*MA,1,(*(*A).A21).VData+j*M21,1);
    }
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

  (*(*A).A12).Type=Admissible;

  (*(*A).A12).A11=NULL;
  (*(*A).A12).A21=NULL;
  (*(*A).A12).A12=NULL;
  (*(*A).A12).A22=NULL;

  (*(*A).A12).Data=NULL;

  (*(*A).A12).Rank=Rank;

  (*(*A).A12).UData=(double complex *)malloc(N12*Rank*sizeof(double complex));
  assert((*(*A).A12).UData != NULL);
  (*(*A).A12).VData=(double complex *)malloc(M12*Rank*sizeof(double complex));
  assert((*(*A).A12).UData != NULL);

  for (j=0; j < Rank; j++)
    {
      cblas_zcopy(N12,(*A).UData+j*NA,1,(*(*A).A12).UData+j*N12,1);
      cblas_zcopy(M12,(*A).VData+j*MA+M11,1,(*(*A).A12).VData+j*M12,1);
    }
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

  (*(*A).A22).Type=Admissible;

  (*(*A).A22).A11=NULL;
  (*(*A).A22).A21=NULL;
  (*(*A).A22).A12=NULL;
  (*(*A).A22).A22=NULL;

  (*(*A).A22).Data=NULL;

  (*(*A).A22).Rank=Rank;

  (*(*A).A22).UData=(double complex *)malloc(N22*Rank*sizeof(double complex));
  assert((*(*A).A22).UData != NULL);
  (*(*A).A22).VData=(double complex *)malloc(M22*Rank*sizeof(double complex));
  assert((*(*A).A22).UData != NULL);

  for (j=0; j < Rank; j++)
    {
      cblas_zcopy(N22,(*A).UData+j*NA+N11,1,(*(*A).A22).UData+j*N22,1);
      cblas_zcopy(M22,(*A).VData+j*MA+M11,1,(*(*A).A22).VData+j*M22,1);
    }
/*
  Free memory
*/
  free((*A).UData); free((*A).VData);

  (*A).Rank=0;
  (*A).UData=NULL;
  (*A).VData=NULL;

  HCheckM(A);
}
