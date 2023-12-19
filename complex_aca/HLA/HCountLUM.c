/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                03.06.2011                                   *
\*****************************************************************************/
/*
  Last change 03.06.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HCountLUM(CBlock *A,long *NumberHLU,long *NumberDLU)
{
  HCheckM(A);
/*
  Local variables
*/
  long NumHLU,NumDLU;
/*
  Check matrix type
*/
  assert(((*A).Type  == Hierarchical) || ((*A).Type == Dense));

  *NumberHLU=0;
  *NumberDLU=0;

  if((*A).Type  == Hierarchical)
    {
      *NumberHLU+=1;

      HCountLUM((*A).A11,&NumHLU,&NumDLU);
      *NumberHLU+=NumHLU;
      *NumberDLU+=NumDLU;

      HCountLUM((*A).A22,&NumHLU,&NumDLU);
      *NumberHLU+=NumHLU;
      *NumberDLU+=NumDLU;
    }
  else if ((*A).Type == Dense)
    {
      *NumberDLU+=1;
    }
}

