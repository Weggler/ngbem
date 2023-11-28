//------------------------------------------------------------------
//                                                                  
//  routine name  -  BlockHMatrixMultZ_fs
//                                                                  
//------------------------------------------------------------------
//
//  last revision -  Aug 12
//
//  purpose       -  matrix vector multiplication of a
//                   MBH-matrix (matrix-driven Block-H-Matrix)
//                   of dimension ndofmaca
//                   resulting from vector-valued fundamental solution
//                   with a vector of dimension ndofmaca
//
//  REMARK: this has nothing to do with the DC-BEM, 
//          i.e., no squeezing or expanding of the
//          rhs and solution vector necessary !
//  
//
//  in:
//      BlockHMat - MBH-Matrix
//              X - vector
//
//  out:          
//              Y == BlockHMat X
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/GlobalVariables.h"

void BlockMHMatrixMultZ_fs(double complex *X,double complex *Y,BlockHMatrix *BlockHMat)
{
  BlockHMatrixMult(*BlockHMat,X,Y);
}
