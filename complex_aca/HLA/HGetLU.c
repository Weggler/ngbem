/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                29.11.2010                                   *
\*****************************************************************************/
/*
  Last change 07.06.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HGetLU(CBlock *A,double HLUEps)
{
  HCheckM(A);
/*
  Local variables
*/
  long NA,MA;
  double complex CM1;
  CBlock *B;
  double AEntry,LAEntry,RCOND,ANORM;
 
  static long CountHLU=0,CountDLU=0;

  static double AIniEntry;
  static long Count=0;
/*
  LAPACK variables
*/
//  MKL_INT NL,*IPIV,INFO,LWORKL;
//  MKL_Complex16 *AL,*WORK;
  int NL,*IPIV,INFO,LWORKL;
  double *RWORK;
  complex double *AL,*WORK;
/*
  Check matrix type
*/
  assert(((*A).Type  == Hierarchical) || ((*A).Type == Dense));
/*
  Dimension
*/
  NA=(*A).NBlock;
  MA=(*A).MBlock;
  assert(NA == MA);

  if (Count == 0)
    {
      AIniEntry=HNormM(A)/sqrt(((double)NA)*((double)MA));
      Count=1;
    }

  if((*A).Type  == Hierarchical)
    {
/*
  Block LU decomposition
*/
      CM1=-1.0+I*0.0;

      HGetLU((*A).A11,HLUEps);

      HSolveLURight((*A).A11,(*A).A21,HLUEps);

      B=(CBlock *)malloc(1*sizeof(CBlock));
      assert(B != NULL);
      HIniM(B);

      HMultMM((*A).A21,(*A).A12,B,HLUEps);

      HScalM(&CM1,B);

      HSumMM(B,(*A).A22,HLUEps);
/*
  Free memory
*/
      HFreeM(B);
      free(B);

      HGetLU((*A).A22,HLUEps);

      HSolveLULeft((*A).A11,(*A).A12,HLUEps);

      CountHLU+=1;
    }
  else if ((*A).Type == Dense)
    {
/*
  LAPACK LU decomposition
*/
    LAEntry=HNormM(A)/sqrt((double)(NA)*(double)(MA));
//      NL=(MKL_INT)NA;
//      AL=(MKL_Complex16 *)(*A).Data;
      NL=(int)NA;
      AL=(complex double *)(*A).Data;

//      IPIV=(MKL_INT *)malloc(NL*sizeof(MKL_INT));
      IPIV=(int *)malloc(NL*sizeof(int));
      assert(IPIV !=NULL);

//      WORK=(MKL_Complex16 *)malloc(10*NL*sizeof(MKL_Complex16));
      WORK=(complex double *)malloc(10*NL*sizeof(complex double));
      assert(WORK != NULL);

      RWORK=(double *)malloc(4*NL*sizeof(double));
      assert(RWORK != NULL);

      LWORKL=10*NL;
 
      zgetrf_(&NL,&NL,AL,&NL,IPIV,&INFO);
      assert(INFO == 0);
/*
  LAPACK inverse
*/
      zgetri_(&NL,AL,&NL,IPIV,WORK,&LWORKL,&INFO);
      assert(INFO == 0);
     
      free(IPIV); free(RWORK); free(WORK);

      CountDLU+=1;
  }

  HCheckM(A);

  if ((*A).Type  == Hierarchical)
    AEntry=HNormM((*A).A12)/sqrt(((double)(*(*A).A12).NBlock)*((double)(*(*A).A12).MBlock));
  else
    AEntry=HNormM(A)/sqrt((double)(NA)*(double)(MA));

  if (AEntry > AIniEntry/(0.01*HLUEps))
    {
      printf("\nElement growth during the HLU : %15.8e instead of %15.8e\n",AEntry,AIniEntry);
      exit(0);
    }

  if (InfoLevel == 1)
    {
      printf("    Done  : Hierarchical LU's %8ld and Direct LU's %8ld\r",CountHLU,CountDLU);
      fflush(stdout);
    }
}
