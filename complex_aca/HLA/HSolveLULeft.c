/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                29.11.2010                                   *
\*****************************************************************************/
/*
  Last change 14.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HSolveLULeft(CBlock *A,CBlock *B,double HLUEps)
{
  HCheckM(A);
  HCheckM(B);
/*
  Local variables
*/
  long NA,MA,NB,MB,N11;
  double complex CM1;
  CBlock *B11,*B21,*C;
/*
  LAPACK variables
*/
  int NL,*IPIV,INFO;
/*
  Check matrix type
*/
  assert(((*A).Type == Hierarchical) || ((*A).Type == Dense));
/*
  Help blocks
*/
  C=(CBlock *)malloc(1*sizeof(CBlock));
  assert(C != NULL);
  HIniM(C);

  B11=(CBlock *)malloc(1*sizeof(CBlock));
  assert(B11 != NULL);
  HIniM(B11);

  B21=(CBlock *)malloc(1*sizeof(CBlock));
  assert(B21 != NULL);
  HIniM(B21);
/*
  Dimension of the matrix A
*/
  NA=(*A).NBlock;
  MA=(*A).MBlock;
  assert(NA == MA);
/*
  Dimension of the matrix B
*/
  NB=(*B).NBlock;
  MB=(*B).MBlock;
  assert(MA == NB);

  CM1=-1.0+I*0.0;

  if (((*A).Type == Hierarchical) && ((*B).Type == Hierarchical))
    {
/*
  L^(-1) Block 2,1
*/
      HMultMM((*A).A21,(*B).A11,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,(*B).A21,HLUEps);
      HFreeM(C);
/*
  L^(-1) Block 2,2
*/
      HMultMM((*A).A21,(*B).A12,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,(*B).A22,HLUEps);
      HFreeM(C);
/*
  D^(-1) Block 1,1
*/
      HSolveLULeft((*A).A11,(*B).A11,HLUEps);
      HSolveLULeft((*A).A11,(*B).A12,HLUEps);
/*
  D^(-1) Block 2,2
*/
      HSolveLULeft((*A).A22,(*B).A21,HLUEps);
      HSolveLULeft((*A).A22,(*B).A22,HLUEps);
/*
  U^(-1) Block 1,1
*/
      HMultMM((*A).A12,(*B).A21,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,(*B).A11,HLUEps);
      HFreeM(C);
/*
  U^(-1) Block 1,2
*/
      HMultMM((*A).A12,(*B).A22,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,(*B).A12,HLUEps);
      HFreeM(C);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Dense))
    {
      N11=(*(*A).A11).MBlock;
      HSplitDenseH(B,N11,B11,B21);
      HFreeM(B);
/*
  L^(-1)
*/
      HMultMM((*A).A21,B11,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,B21,HLUEps);
      HFreeM(C);
/*
  D^(-1)
*/
      HSolveLULeft((*A).A11,B11,HLUEps);
      HSolveLULeft((*A).A22,B21,HLUEps);
/*
  U^(-1)
*/
      HMultMM((*A).A12,B21,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,B11,HLUEps);
      HFreeM(C);

      HJoinDenseH(B,B11,B21,HLUEps);

      HFreeM(B11); HFreeM(B21);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Admissible))
    {
      (*C).IBlock=(*A).IBlock;
      (*C).JBlock=(*B).JBlock;
      (*C).NBlock=NB;
      (*C).MBlock=(*B).Rank;

      (*C).Type=Dense;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=(double complex *)malloc(NB*((*B).Rank)*sizeof(double complex));
      assert((*C).Data != NULL);

      cblas_zcopy(NB*((*B).Rank),(*B).UData,1,(*C).Data,1);

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;

      HSolveLULeft(A,C,HLUEps);

      if ((*C).Type == Admissible)
	HAToD(C);

      cblas_zcopy(NB*((*B).Rank),(*C).Data,1,(*B).UData,1);

      HFreeM(C);
    }
  else
    {
      HMultMM(A,B,C,HLUEps);
      HFreeM(B);
      HCopyM(C,B);
      HFreeM(C);
    }
/*
  Free memory
*/
  free(C); free(B11); free(B21);

  HCheckM(A);
  HCheckM(B);   
}
