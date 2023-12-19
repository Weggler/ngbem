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

void HJoinDenseV(CBlock *A,CBlock *A11,CBlock *A12,double HLUEps)
{

  HCheckM(A);
  HCheckM(A11);
  HCheckM(A12);
/*
  Local variables
*/
  long NA,MA,N11,M11,N12,M12;

  HCheckM(A);
  HCheckM(A11);
  HCheckM(A12);

  if ((*A11).Type == Null)
    HNToD(A11);
  else if ((*A11).Type == Hierarchical)
    HHToD(A11);
  else if ((*A11).Type == Admissible)
    HAToD(A11);

  if ((*A12).Type == Null)
    HNToD(A12);
  else if ((*A12).Type == Hierarchical)
    HHToD(A12);
  else if ((*A12).Type == Admissible)
    HAToD(A12);

  N11=(*A11).NBlock;
  M11=(*A11).MBlock;

  N12=(*A12).NBlock;
  M12=(*A12).MBlock;

  assert(N11 == N12);
  assert((M11 >= 1) && (M12 >=1));

  NA=N11;
  MA=M11+M12;

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

  cblas_zcopy(N11*M11,(*A11).Data,1,(*A).Data,1);
  cblas_zcopy(N12*M12,(*A12).Data,1,(*A).Data+N11*M11,1);

  (*A).Rank=0;
  (*A).UData=NULL;
  (*A).VData=NULL;

  HDToA(A,HLUEps);

  HCheckM(A);
  HCheckM(A11);
  HCheckM(A12);
}
