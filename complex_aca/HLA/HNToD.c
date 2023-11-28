/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                08.12.2010                                   *
\*****************************************************************************/
/*
  Last change 08.12.2010
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HNToD(CBlock *A)
{
/*
  Local variables
*/
  long NA,MA;
  double complex C0,C1;

  assert((*A).Type == Null);

  (*A).Type=Dense;

  NA=(*A).NBlock;
  MA=(*A).MBlock;

  (*A).Data=(double complex *)malloc(NA*MA*sizeof(double complex));
  assert((*A).Data != NULL);

  zset(NA*MA,0.0+I*0.0,(*A).Data,1);
}
