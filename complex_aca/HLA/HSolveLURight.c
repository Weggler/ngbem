/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                01.12.2010                                   *
\*****************************************************************************/
/*
  Last change 14.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HSolveLURight(CBlock *A,CBlock *B,double HLUEps)
{
  HCheckM(A);
  HCheckM(B);
/*
  Local variables
*/
  long NA,MA,NB,MB,M11;
  long i,j;
  double complex CM1;
  CBlock *B11,*B12,*C;
/*
  LAPACK variables
*/
//  MKL_Complex16 *AL;
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

  B12=(CBlock *)malloc(1*sizeof(CBlock));
  assert(B12 != NULL);
  HIniM(B12);
/*
  Dimension of the matrix A
*/
  NA=(*A).NBlock;
  MA=(*A).MBlock;
//  printf("HSolveLURight: NA = %d, MA = %d\n", NA, MA);
  assert(NA == MA);
/*
  Dimension of the matrix B
*/
  NB=(*B).NBlock;
  MB=(*B).MBlock;
//  printf("HSolveLURight: NA = %d, NB = %d, MB = %d\n", NA, NB, MB);
  assert(NA == MB);

  CM1=-1.0+I*0.0;

  if (((*A).Type == Hierarchical) && ((*B).Type == Hierarchical))
    {
/*
  U^(-1) Block 1,2
*/
      HMultMM((*B).A11,(*A).A12,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,(*B).A12,HLUEps);
      HFreeM(C);
/*
  U^(-1) Block 2,2
*/
      HMultMM((*B).A21,(*A).A12,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,(*B).A22,HLUEps);
      HFreeM(C);
/*
  D^(-1) Block 1,1
*/
     HSolveLURight((*A).A11,(*B).A11,HLUEps);
     HSolveLURight((*A).A11,(*B).A21,HLUEps);
/*
  D^(-1) Block 2,2
*/
      HSolveLURight((*A).A22,(*B).A12,HLUEps);
      HSolveLURight((*A).A22,(*B).A22,HLUEps);
/*
  L^(-1) Block 1,1
*/
      HMultMM((*B).A12,(*A).A21,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,(*B).A11,HLUEps);
      HFreeM(C);
/*
  L^(-1) Block 2,1
*/
      HMultMM((*B).A22,(*A).A21,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,(*B).A21,HLUEps);
      HFreeM(C);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Dense))
    {
      M11=(*(*A).A11).NBlock;
      HSplitDenseV(B,M11,B11,B12);
      HFreeM(B);
/*
  U^(-1)
*/
      HMultMM(B11,(*A).A12,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,B12,HLUEps);
      HFreeM(C);
/*
  D^(-1)
*/
      HSolveLURight((*A).A11,B11,HLUEps);
      HSolveLURight((*A).A22,B12,HLUEps);
/*
  L^(-1)
*/
      HMultMM(B12,(*A).A21,C,HLUEps);
      HScalM(&CM1,C);
      HSumMM(C,B11,HLUEps);
      HFreeM(C);
 
      HJoinDenseV(B,B11,B12,HLUEps);

      HFreeM(B11); HFreeM(B12);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Admissible))
    {
      (*C).IBlock=(*B).IBlock;
      (*C).JBlock=(*A).JBlock;
      (*C).NBlock=(*B).Rank;
      (*C).MBlock=MB;

      (*C).Type=Dense;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=(double complex *)malloc(((*B).Rank)*MB*sizeof(double complex));
      assert((*C).Data != NULL);

      for (j=0; j < MB; j++)
       for (i=0; i < (*B).Rank ; i++)
	 (*C).Data[j*(*B).Rank+i]=conj((*B).VData[i*MB+j]);

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;

      HSolveLURight(A,C,HLUEps);

      if ((*C).Type == Admissible)
	HAToD(C);

      for (j=0; j <  (*B).Rank; j++)
       for (i=0; i < MB ; i++)
	 (*B).VData[j*MB+i]=conj((*C).Data[i*(*B).Rank+j]);

      HFreeM(C);
    }
  else
    {
      HMultMM(B,A,C,HLUEps);
      HFreeM(B);
      HCopyM(C,B);
      HFreeM(C);
    }
/*
  Free memory
*/
  HFreeM(C); HFreeM(B11); HFreeM(B12);

  free(C); free(B11); free(B12);

  HCheckM(A);
  HCheckM(B);   
}
