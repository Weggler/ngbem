/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                26.01.2011                                   *
\*****************************************************************************/
/*
  Last change 26.01.2011
*/
/*
  Includes for all functions 
*/

# include "Includes.h"

void HHHHBlock(long IB,long JB, long N, long M,CBlock *P,CBlock *A)
{
  A=(CBlock *)malloc(1*sizeof(CBlock));
  assert(A != NULL);
/*
  Preparations
*/
  (*A).IBlock=IB;
  (*A).JBlock=JB;
  (*A).NBlock=N;
  (*A).MBlock=M;

  (*A).Type=Hierarchical;

  (*A).Data=NULL;

  (*A).Rank=0;
  (*A).UData=NULL;
  (*A).VData=NULL;
/*
  Blocks
*/
  (*A).A11=P[0];
  (*A).A21=P[1];
  (*A).A12=P[2];
  (*A).A22=P[3];
}
