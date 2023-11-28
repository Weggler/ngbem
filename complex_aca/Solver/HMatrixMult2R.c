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

# include "Includes.h"

void HMatrixMult2R(double *X,double *Y,HMatrix *HMat)
{
/*
  Local variables
*/
  long N,M,i,j,IBlock,JBlock,NBlock,MBlock,Rank;
  double complex *XCopy,*YCopy,*Work,C0,C1;
  BlockType BlType;
/*
  Initialisation
*/
  C0=0.0+0.0*I;
  C1=1.0+0.0*I;

  N=HMat->NRow; M=HMat->NColumn;
/*
  Memory allocation for the copies
*/
  XCopy=(double  complex *)malloc(M*sizeof(double complex));
  assert(XCopy != NULL);

  YCopy=(double  complex *)malloc(N*sizeof(double complex));
  assert(YCopy != NULL);
/*
  Permutation of columns
*/
  for (j=0; j < M; j++)
    XCopy[j]=X[HMat->PermuColumn[j]]+X[M+HMat->PermuColumn[j]]*I;

  zset(N,0.0+0.0*I,YCopy,1);
/*
  All blocks
*/
  for (i=0; i < HMat->NBlocks; i++)
    {
      BlType=HMat->CBlocks[i].Type;

      IBlock=HMat->CBlocks[i].IBlock;
      JBlock=HMat->CBlocks[i].JBlock;

      NBlock=HMat->CBlocks[i].NBlock;
      MBlock=HMat->CBlocks[i].MBlock;

      if (BlType == Dense)
	{
	  cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,MBlock,
		      &C1,HMat->CBlocks[i].Data,NBlock,
		      XCopy+JBlock,1,&C1,YCopy+IBlock,1);

	  if(HMat->Symmetry == ComplexSymmetric)
	    cblas_zgemv(CblasColMajor,CblasTrans,NBlock,MBlock,
			&C1,HMat->CBlocks[i].Data,NBlock,
			XCopy+IBlock,1,&C1,YCopy+JBlock,1);
	}
      else if (BlType == Admissible)
	{
	  Rank=HMat->CBlocks[i].Rank;

	  Work=(double complex *)malloc(Rank*sizeof(double complex));
	  assert(Work != NULL);

	  cblas_zgemv(CblasColMajor,CblasConjTrans,MBlock,Rank,
		      &C1,HMat->CBlocks[i].VData,MBlock,
		      XCopy+JBlock,1,&C0,Work,1);

	  cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,Rank,
		      &C1,HMat->CBlocks[i].UData,NBlock,
		      Work,1,&C1,YCopy+IBlock,1);

	  if(HMat->Symmetry == ComplexSymmetric)
	    {
	      cblas_zgemv(CblasColMajor,CblasTrans,NBlock,Rank,
			  &C1,HMat->CBlocks[i].UData,NBlock,
			  XCopy+JBlock,1,&C0,Work,1);

	      cblas_zgemv(CblasRowMajor,CblasConjTrans,MBlock,Rank,
			  &C1,HMat->CBlocks[i].VData,Rank,
			  Work,1,&C1,YCopy+IBlock,1);
	    }

          free(Work); 
	}
    }

  for (i=0; i < N; i++)
    {
      Y[HMat->PermuRow[i]]=creal(YCopy[i]);
      Y[N+HMat->PermuRow[i]]=cimag(YCopy[i]);
    }

  free(XCopy); free(YCopy); 
}
