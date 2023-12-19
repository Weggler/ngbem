/*------------------------------------------------------------------*/
/*                                                                  */
/* S. Rjasanow: Adaptive Cross Approximation                        */
/*                                                                  */
/*  routine name  - HLUSolveMM                                      */
/*                                                                  */
/*------------------------------------------------------------------*/
/*                                                                  */
/*  last revision -  Feb 11                                         */
/*  purpose       -  LU solver for hierarchic matrix                */
/*                                                                  */
/*  in                                                              */
/*           NRHS - number of right hand sides                      */
/*           HMat - hierarchic matrix                               */
/*              Y - right hand side                                 */
/*         HLUEps - accuracy of LU decomposition                    */
/*                                                                  */
/*  out                                                             */
/*              X - solution vector                                 */
/*                                                                  */
/*------------------------------------------------------------------*/

/* Includes for all functions */
# include "../CInclude/Includes.h"

void HLUSolveMM(long NRHS,HMatrix *HMat,
                double complex *X,
                double complex *Y,
                double HLUEps)
{

  /* Local variables */
  long i,j,N,M;
  CBlock *BY;

  /* Dimension of the matrix */
  N=(*HMat).NRow; 
  M=(*HMat).NColumn;

  assert(N == M);

  /* Memory allocation for the RHS */
  BY=(CBlock *)malloc(1*sizeof(CBlock));
  assert(BY != NULL);
  HIniM(BY);

  /* Preparations */
  (*BY).IBlock=1;
  (*BY).JBlock=1;
  (*BY).NBlock=N;
  (*BY).MBlock=NRHS;

  (*BY).Type=Dense;

  (*BY).A11=NULL;
  (*BY).A21=NULL;
  (*BY).A12=NULL;
  (*BY).A22=NULL;

  (*BY).Data=(double complex *)malloc(N*NRHS*sizeof(double complex));
  assert((*BY).Data != NULL);

  /* Permutation of the RHS */
  for (j=0; j < NRHS; j++)
  {
    for (i=0; i < N; i++)
    {
      (*BY).Data[j*N+i]=Y[j*N+(*HMat).PermuRow[i]];
    }
  }

  (*BY).Rank=0;
  (*BY).UData=NULL;
  (*BY).VData=NULL;

  /* LU */
  HGetLU((*HMat).CBlocks,HLUEps);

  /* Solution */
  HSolveLULeft((*HMat).CBlocks,BY,HLUEps);

  /* Permutation of of the solution */
  for (j=0; j < NRHS; j++)
    for (i=0; i < N; i++)
 	X[j*N+(*HMat).PermuColumn[i]]=(*BY).Data[j*N+i];

  /* Free memory */
  HFreeM(BY);

  free(BY);
}
