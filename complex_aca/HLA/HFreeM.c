/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                28.11.2010                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HFreeM(CBlock *A)
{
  HCheckM(A);

  if((*A).Type == Hierarchical)
    {
      HFreeM((*A).A11);
      HFreeM((*A).A21);
      HFreeM((*A).A12);
      HFreeM((*A).A22);

      free((*A).A11);
      (*A).A11=NULL;

      free((*A).A21);
      (*A).A21=NULL;

      free((*A).A12);
      (*A).A12=NULL;

      free((*A).A22);
      (*A).A22=NULL;
    }
  else if((*A).Type == Dense)
    {
      free((*A).Data);

      (*A).Data=NULL;
    }
  else if ((*A).Type == Admissible)
    {
      free((*A).UData);
      free((*A).VData);

      (*A).Rank=0;
      (*A).UData=NULL;
      (*A).VData=NULL;
    }
   
  HIniM(A);

  HCheckM(A);
}
