/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                  Adaptive Cross Approximation for EMSS                      *
*                                                                             *
*                                26.01.2011                                   *
\*****************************************************************************/
/*
  Last change 04.06.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/GlobalVariables.h"

void BlockHLUSolveMM(long NRHS, 
                     BlockHMatrix *BHMat, 
                     double complex *X, 
                     double complex *Y, 
                     double HLUEps)
{
/*
  Local variables
*/
  long NB,MB,N,M,NumHLU,NumDLU;
  double MemA;
  CBlock *A,*BY;
/*
  Memory allocation
*/
  BY=(CBlock *)malloc(1*sizeof(CBlock));
  assert(BY != NULL);
/*
  Transform the matrix in suitable format of blocks of 2 by 2 matrices, 
  where filled up with zero and unity submatrices if necessary...
*/
  BHMatrixToBlock(BHMat,&A);
/*
  Transform the right hand side, same here
*/
  RHSToBlock(NRHS,BHMat,Y,BY);
/*
  LU
*/
  if (InfoLevel == 1)
    {
     MemA=HMemM(A);
     HCountLUM(A,&NumHLU,&NumDLU);

      printf("  Start hierarchical LU decomposition with %8.2f MB\n",MemA/(1024.0*1024.0));
      printf("    To do : Hierarchical LU's %8ld and Direct LU's %8ld\n",NumHLU,NumDLU);
      fflush(stdout);
    }

  HGetLU(A,HLUEps);

  if (InfoLevel == 1)
    {
      MemA=HMemM(A);

      printf("    Done with %8.2f MB\n",MemA/(1024.0*1024.0));
      fflush(stdout);
    }
/*
  Solution
*/
  HSolveLULeft(A,BY,HLUEps);
/*
  Transform the solution
*/
  BlockToSolution(NRHS,BHMat,X,BY);
/*
  Free memory
*/
//  HFreeM(A);
  free(A); 

  HFreeM(BY);
  free(BY);
}
