/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                11.11.2010                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HScalM(double complex *Alpha,CBlock *A)
{
  HCheckM(A);

  if((*A).Type  == Hierarchical)
    {
      HScalM(Alpha,(*A).A11);
      HScalM(Alpha,(*A).A21);
      HScalM(Alpha,(*A).A12);
      HScalM(Alpha,(*A).A22);
    }
  else if ((*A).Type == Dense)
    cblas_zscal(((*A).NBlock)*((*A).MBlock),Alpha,(*A).Data,1);
  else if ((*A).Type == Admissible)
    cblas_zscal(((*A).NBlock)*((*A).Rank),Alpha,(*A).UData,1);

  HCheckM(A);
}
