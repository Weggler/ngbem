/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                29.01.2011                                   *
\*****************************************************************************/
/*
  Last change 29.01.2011
*/
/*
  Includes for all functions 
*/

# include "Includes.h"

void HNNN4x4Block(long IB,long JB, long N, long M,CBlock *P,CBlock *A)
{
/*
  Local variables
*/
  long I11,I21,I12,I22,J11,J21,J12,J22, N11,N21,N12,N22,M11,M21,M12,M22;
  CBlock PP[4];
 
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
  (*A).A11=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A11 != NULL);
  (*A).A21=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A21 != NULL);
  (*A).A12=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A12 != NULL);
  (*A).A22=(CBlock *)malloc(1*sizeof(CBlock));
  assert((*A).A22 != NULL);
/*
  Positions
*/
  I11=IB;
  J11=JB;

  I21=I11+2;
  J21=J11;

  I12=I11;
  J12=J11+P[0].MBlock+P[1].MBlock;

  I22=I11+2;
  J22=J11+P[0].MBlock+P[1].MBlock;
/*
  Dimesions
*/
  N11=P[0].NBlock+1;
  M11=P[0].MBlock+P[1].MBlock;

  N21=2;
  M21=M11;

  N12=P[2].MBlock+1;
  M12=P[2].MBlock+P[3].MBlock;

  assert(N12 == N11);

  N22=2;
  M22=M12;
/*
  Block 1,1
*/
  PP[0]=P[0];
  PP[1]=NULL;
  PP[2]=P[1];
  PP[3]=NULL;

  HNHNBlock(I11,J11,N11,M11,PP,(*A).A11);
/*
  Block 2,1
*/
  HNBlock(I21,J21,N21,M21,PP,(*A).A21);
/*
  Block 1,2
*/
  PP[0]=P[2];
  PP[1]=NULL;
  PP[2]=P[3];
  PP[3]=NULL;

  HNHNBlock(I12,J12,N12,M12,PP,(*A).A12);
/*
  Block 2,2
*/
  HNBlock(I22,J22,N22,M22,PP,(*A).A22);
}
