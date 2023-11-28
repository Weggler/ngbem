/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                16.03.2010                                   *
\*****************************************************************************/
/*
  Last change 16.03.2010
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void GMResPrint(int Iter,double Error)
{
  printf("  Iter = %8d, Error = %15.8le\r",Iter,Error);
  fflush(stdout);
}
