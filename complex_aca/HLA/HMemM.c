/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                19.01.2011                                   *
\*****************************************************************************/
/*
  Last change 19.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

double HMemM(CBlock *A)
{
  double Memory;

  HCheckM(A);

  Memory=0.0;

  if ((*A).Type  == Null)
    Memory+=(double)sizeof(CBlock);
  else if ((*A).Type  == Hierarchical)
    {
      Memory+=(HMemM((*A).A11)+HMemM((*A).A21)+HMemM((*A).A12)+HMemM((*A).A22));
      Memory+=4.0*(double)sizeof(CBlock);
    }
  else if ((*A).Type == Dense)
    Memory+=(double)((*A).NBlock*(*A).MBlock*sizeof(double complex));
  else if ((*A).Type == Admissible)
    Memory+=(double)((*A).Rank*((*A).NBlock+(*A).MBlock)*sizeof(double complex));

  return Memory;
}
