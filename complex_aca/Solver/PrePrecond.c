/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                08.05.2010                                   *
\*****************************************************************************/
/*
  Last change 08.05.2010
*/
/*
  Includes for all functions 
*/

# include "Includes.h"

void PrePrecond(HMatrix *HMat,double complex (*ElementName)(long,long),double complex *y)
{
/*
  Local variables
*/
  long i,j,N,Error;
  PrecondType Precond;
  double complex *XCopy,*YCopy;

  N=HMat->NRow;
  Precond=HMat->Precond;

  switch(Precond)
    {
    case None:
      printf("In None\n");

      break;
    case Jacobi:

      HMat->PrecondIData=NULL;
      HMat->PrecondDData=(double complex *)malloc(N*sizeof(double complex));

      for (i=0; i < N; i++)
	{
	  (HMat->PrecondDData)[i]=ElementName(i,i);

	  y[i]=y[i]/(HMat->PrecondDData)[i];
	}
      break;
    case TriDiag:
      printf("In TriDiag\n");

      HMat->PrecondIData=NULL;
      HMat->PrecondDData=(double complex *)malloc(3*N*sizeof(double complex));
      assert(HMat->PrecondDData != NULL);

      (HMat->PrecondDData)[0]=0.0+0.0*I;
      (HMat->PrecondDData)[N+0]=ElementName(0,1);
      (HMat->PrecondDData)[2*N+0]=ElementName(0,0);

      for (i=1; i < N-1; i++)
	{
	  (HMat->PrecondDData)[i]=ElementName(i-1,i);
	  (HMat->PrecondDData)[N+i]=ElementName(i,i+1);
	  (HMat->PrecondDData)[2*N+i]=ElementName(i,i);
	}

      (HMat->PrecondDData)[N-1]=ElementName(N-2,N-1);
      (HMat->PrecondDData)[2*N-1]=0.0+0.0*I;
      (HMat->PrecondDData)[3*N-1]=ElementName(N-1,N-1);

      SolveTriDiag(N-1,(HMat->PrecondDData),(HMat->PrecondDData)+N,(HMat->PrecondDData)+2*N,y,y);       
      break;
    case ILU:
      printf("In ILU\n");

      HMat->PrecondIData=NULL;
      HMat->PrecondDData=NULL;

      GetSparseStructure(HMat);

      XCopy=(double complex *)malloc(N*sizeof(double complex));
      assert(XCopy != NULL);

      YCopy=(double complex *)malloc(N*sizeof(double complex));
      assert(YCopy != NULL);

      for (j=0; j < N; j++)
	YCopy[j]=y[(HMat->PermuColumn)[j]];

      Error=zlusolC(YCopy,XCopy,HMat->lu);
      assert(Error == 0);

      for (i=0; i < N; i++)
	y[(HMat->PermuRow)[i]]=XCopy[i];

      free(XCopy); free(YCopy);
      break;
    }
}
