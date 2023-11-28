/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                08.12.2010                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HHToD(CBlock *A)
{
  HCheckM(A);
/*
  Local variables
*/
  long NA,MA,N11,M11,N21,M21,N12,M12,N22,M22;
  long i,j;

  assert((*A).Type == Hierarchical);

  (*A).Type=Dense;

  NA=(*A).NBlock;
  MA=(*A).MBlock;

  (*A).Data=(double complex *)malloc(NA*MA*sizeof(double complex));
  assert((*A).Data != NULL);
/*
  Block 1,1
*/
  N11=(*(*A).A11).NBlock;
  M11=(*(*A).A11).MBlock;

  if ((*(*A).A11).Type == Hierarchical)
    HHToD((*A).A11);
  else if ((*(*A).A11).Type == Admissible)
    HAToD((*A).A11);
  else if ((*(*A).A11).Type == Null)
    HNToD((*A).A11);

  for (j=0; j < M11; j++)
    cblas_zcopy(N11,(*(*A).A11).Data+j*N11,1,(*A).Data+j*NA,1);
/*
  Block 2,1
*/
  N21=NA-N11;
  M21=M11;

  if ((*(*A).A21).Type == Hierarchical)
    HHToD((*A).A21);
  else if ((*(*A).A21).Type == Admissible)
    HAToD((*A).A21);
  else if ((*(*A).A21).Type == Null)
    HNToD((*A).A21);

  for (j=0; j < M21; j++)
    cblas_zcopy(N21,(*(*A).A21).Data+j*N21,1,(*A).Data+j*NA+N11,1);
/*
  Block 1,2
*/
  N12=N11;
  M12=MA-M11;

  if ((*(*A).A12).Type == Hierarchical)
    HHToD((*A).A12);
  else if ((*(*A).A12).Type == Admissible)
    HAToD((*A).A12);
  else if ((*(*A).A12).Type == Null)
    HNToD((*A).A12);

  for (j=0; j < M12; j++)
    cblas_zcopy(N12,(*(*A).A12).Data+j*N12,1,(*A).Data+(j+M11)*NA,1);
/*
  Block 2,2
*/
  N22=NA-N11;
  M22=MA-M11;
 
  if ((*(*A).A22).Type == Hierarchical)
    HHToD((*A).A22);
  else if ((*(*A).A22).Type == Admissible)
    HAToD((*A).A22);
  else if ((*(*A).A22).Type == Null)
    HNToD((*A).A22);

  for (j=0; j < M22; j++)
    cblas_zcopy(N22,(*(*A).A22).Data+j*N22,1,(*A).Data+(j+M11)*NA+N11,1);
/*
  Free blocks
*/
  HFreeM((*A).A11);
  HFreeM((*A).A21);
  HFreeM((*A).A12);
  HFreeM((*A).A22);

  free((*A).A11); 
  free((*A).A21); 
  free((*A).A12); 
  free((*A).A22);

  (*A).A11=NULL;
  (*A).A21=NULL;
  (*A).A12=NULL;
  (*A).A22=NULL;

  HCheckM(A);
}
