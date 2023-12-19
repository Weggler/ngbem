/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                10.12.2010                                   *
\*****************************************************************************/
/*
  Last change 18.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HJoinDenseH(CBlock *A,CBlock *A11,CBlock *A21,double HLUEps)
{

  HCheckM(A);
  HCheckM(A11);
  HCheckM(A21);
/*
  Local variables
*/
  long NA,MA,N11,M11,N21,M21;
  long j;

  if ((*A11).Type == Null)
    HNToD(A11);
  else if ((*A11).Type == Hierarchical)
    HHToD(A11);
  else if ((*A11).Type == Admissible)
    HAToD(A11);

  if ((*A21).Type == Null)
    HNToD(A21);
  else if ((*A21).Type == Hierarchical)
    HHToD(A21);
  else if ((*A21).Type == Admissible)
    HAToD(A21);


  N11=(*A11).NBlock;
  M11=(*A11).MBlock;

  N21=(*A21).NBlock;
  M21=(*A21).MBlock;

  assert(M11 == M21);
  assert((N11 >= 1) && (N21 >=1));

  NA=N11+N21;
  MA=M11;

  (*A).IBlock=(*A11).IBlock;
  (*A).JBlock=(*A11).JBlock;
  (*A).NBlock=NA;
  (*A).MBlock=MA;

  (*A).Type=Dense;

  (*A).A11=NULL;
  (*A).A21=NULL;
  (*A).A12=NULL;
  (*A).A22=NULL;

  (*A).Data=(double complex *)malloc(NA*MA*sizeof(double complex));
  assert((*A).Data != NULL);

  for (j=0; j < MA; j++)
    {
      cblas_zcopy(N11,(*A11).Data+j*N11,1,(*A).Data+j*NA,1);
      cblas_zcopy(N21,(*A21).Data+j*N21,1,(*A).Data+j*NA+N11,1);
    }

  HDToA(A,HLUEps);

  HCheckM(A);
  HCheckM(A11);
  HCheckM(A21);
}
