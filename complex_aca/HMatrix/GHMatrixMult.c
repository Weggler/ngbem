/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                29.04.2009                                   *
\*****************************************************************************/
/*
  Last change 29.04.2009
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void GHMatrixMult(GHMatrix GHMat,double complex *X,double complex *Y)
{
/*
  Local variables
*/
  long N,M,i,j,NBlock,MBlock,Rank;
  double complex *XCopy,*YCopy,*Work,C0,C1;
  BlockType BlType;
/*
  Initialisation
*/
  C0=0.0+0.0*I;
  C1=1.0+0.0*I;

  N=GHMat.NRow; M=GHMat.NColumn;
/*
  Memory allocation for the copies
*/
  XCopy=(double complex *)malloc(GHMat.MMax*sizeof(double complex));
  assert(XCopy != NULL);

  YCopy=(double complex *)malloc(GHMat.NMax*sizeof(double complex));
  assert(YCopy != NULL);

  zset(N,C0,Y,1);
/*
  All blocks
*/
  for (i=0; i < GHMat.NBlocks; i++)
    {
      BlType=GHMat.GCBlocks[i].Type;

      NBlock=GHMat.GCBlocks[i].NBlock;
      MBlock=GHMat.GCBlocks[i].MBlock;

/*
  Copy X
*/
      for (j=0; j < MBlock; j++)
	XCopy[j]=X[GHMat.GCBlocks[i].Cols[j]];

      if (BlType == Dense)
	cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,MBlock,
		    &C1,GHMat.GCBlocks[i].Data,NBlock,
		    XCopy,1,&C0,YCopy,1);
      else if (BlType == Admissible)
	{
	  Rank=GHMat.GCBlocks[i].Rank;

	  Work=(double complex *)malloc(Rank*sizeof(double complex));
	  assert(Work != NULL);

	  cblas_zgemv(CblasColMajor,CblasConjTrans,MBlock,Rank,
		      &C1,GHMat.GCBlocks[i].VData,MBlock,
		      XCopy,1,&C0,Work,1);

	  cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,Rank,
		      &C1,GHMat.GCBlocks[i].UData,NBlock,
		      Work,1,&C0,YCopy,1);

          free(Work); 
	}
/*
  Add to Y
*/
      for (j=0; j < NBlock; j++)
	Y[GHMat.GCBlocks[i].Rows[j]]+=YCopy[j];
    }

  free(XCopy); free(YCopy); 
}
