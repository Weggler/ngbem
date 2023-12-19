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

void HMatrixMultZ(double complex *X,double complex *Y,HMatrix *HMat)
{
/*
  Local variables
*/
  long i,j,N,Error;
  PrecondType Precond;
  double complex *XCopy,*YCopy;

  HMatrixMult(*HMat,X,Y);

  N=HMat->NRow;
  Precond=HMat->Precond;

  switch(Precond)
    {
    case None:

      break;
    case Jacobi:

      for (i=0; i < N; i++)
	Y[i]=Y[i]/(HMat->PrecondDData)[i];

      break;
    case TriDiag:

      SolveTriDiag(N-1,(HMat->PrecondDData),(HMat->PrecondDData)+N,(HMat->PrecondDData)+2*N,Y,Y);       
      break;
    case ILU:

      XCopy=(double complex *)malloc(N*sizeof(double complex));
      assert(XCopy != NULL);

      YCopy=(double complex *)malloc(N*sizeof(double complex));
      assert(YCopy != NULL);

      for (j=0; j < N; j++)
	YCopy[j]=Y[(HMat->PermuColumn)[j]];

      Error=zlusolC(YCopy,XCopy,HMat->lu);
      assert(Error == 0);

      for (i=0; i < N; i++)
	Y[(HMat->PermuRow)[i]]=XCopy[i];

      free(XCopy); free(YCopy);
      break;
    }
}
