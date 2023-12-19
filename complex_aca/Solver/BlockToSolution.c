/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                16.02.2011                                   *
\*****************************************************************************/
/*
  Last change 16.02.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void BlockToSolution(long NRHS,
                     BlockHMatrix *BHMat,
                     double complex *X,
                     CBlock *BY)
{
/*
  Local variables
*/
  long N,M,NB,MB;
  long N2,ND;
  long jr,jb,i,ia,ie;
/*
  Block dimension of the matrix
*/
  NB=(*BHMat).NBlockRow; 
  MB=(*BHMat).NBlockColumn;
  assert(NB == MB);
/*
  Dimension of the matrix
*/
  N=(*BHMat).NRow; 
  M=(*BHMat).NColumn;
  assert(N == M);
/*
  Next power of 2
*/
  N2=(long)pow(2,(long)(log(NB-0.5)/log(2.0)+1));

  ND=N2-NB;
/*
  Permutation of the solution
*/
  for (jr=0; jr < NRHS; jr++)
  {
    for (jb=0; jb < MB; jb++)
    {
      ia=(*BHMat).NDim[jb];
      if(jb == MB-1)
	    ie=M;
      else
	    ie=(*BHMat).NDim[jb+1];
	  for (i=ia; i < ie; i++)
      {
//        printf("BlockToSolution ia, ie: %d %d\n", ia,ie);
        X[jr*N+ia+((*BHMat).HBlocks[jb*MB]).PermuColumn[i-ia]]=(*BY).Data[jr*(N+ND)+i];
      }
    }
  }
}
