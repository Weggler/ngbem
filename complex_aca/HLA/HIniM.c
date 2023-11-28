/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                07.01.2011                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HIniM(CBlock *A)
{
  (*A).IBlock=0;
  (*A).JBlock=0;
  (*A).NBlock=1;
  (*A).MBlock=1;

  (*A).Type=Null;

  (*A).A11=NULL;
  (*A).A21=NULL;
  (*A).A12=NULL;
  (*A).A22=NULL;

  (*A).Data=NULL;

  (*A).Rank=0;

  (*A).UData=NULL;
  (*A).VData=NULL;
}
