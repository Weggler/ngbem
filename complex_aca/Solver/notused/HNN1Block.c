/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                26.01.2011                                   *
\*****************************************************************************/
/*
  Last change 29.01.2011
*/
/*
  Includes for all functions 
*/

# include "Includes.h"

void HNN1Block(long IB,long JB, long N, long M,CBlock *P,CBlock *A)
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
  Block 1,1
*/
  (*A).A11=P[0];
/*
  Block 2,1
*/
  (*A).A21=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A21 != NULL);

  (*(*A).A21).IBlock=IB+P[0].NBlock;
  (*(*A).A21).JBlock=JB;
  (*(*A).A21).NBlock=1;
  (*(*A).A21).MBlock=P[0].MBlock;

  (*(*A).A21).Type=Null;

  (*(*A).A21).A11=NULL;
  (*(*A).A21).A21=NULL;
  (*(*A).A21).A12=NULL;
  (*(*A).A21).A22=NULL;

  (*(*A).A21).Data=NULL;

  (*(*A).A21).Rank=0;
  (*(*A).A21).UData=NULL;
  (*(*A).A21).VData=NULL;
/*
  Block 1,2
*/
  (*A).A12=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A12 != NULL);

  (*(*A).A12).IBlock=IB;
  (*(*A).A12).JBlock=JB+P[0].MBlock;
  (*(*A).A12).NBlock=P[0].NBlock;
  (*(*A).A12).MBlock=1;

  (*(*A).A12).Type=Null;

  (*(*A).A12).A11=NULL;
  (*(*A).A12).A21=NULL;
  (*(*A).A12).A12=NULL;
  (*(*A).A12).A22=NULL;

  (*(*A).A12).Data=NULL;

  (*(*A).A12).Rank=0;
  (*(*A).A12).UData=NULL;
  (*(*A).A12).VData=NULL;
/*
  Block 2,2
*/
  (*A).A22=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A22 != NULL);

  (*(*A).A22).IBlock=IB+P[0].NBlock;
  (*(*A).A22).JBlock=JB+P[0].MBlock;
  (*(*A).A22).NBlock=1;
  (*(*A).A22).MBlock=1;

  (*(*A).A22).Type=Dense;

  (*(*A).A22).A11=NULL;
  (*(*A).A22).A21=NULL;
  (*(*A).A22).A12=NULL;
  (*(*A).A22).A22=NULL;

  (*(*A).A22).Data=(double complex *)malloc(1*sizeof(double complex));
  assert(*(*A).A22).Data != NULL);

  (*(*A).A22).Data[0]=1.0+I*0.0;

  (*(*A).A22).Rank=0;
  (*(*A).A22).UData=NULL;
  (*(*A).A22).VData=NULL;
}
