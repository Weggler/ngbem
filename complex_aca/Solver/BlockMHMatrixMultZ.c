//------------------------------------------------------------------
//                                                                  
//  routine name  -  BlockMHMatrixMultZ
//                                                                  
//------------------------------------------------------------------
//
//  last revision -  Feb 13
//
//  purpose       -  matrix vector multiplication of a
//                   MBH-matrix (matrix-driven Block-H-Matrix)
//                   with a vector
//                   NOTE: this routine copes with the interplay 
//                         between the discontinous Galerkin
//                         scheme MACA and the conventional, 
//                         continuous Galerkin scheme !
//
//  in:
//      BlockHMat - MBH-Matrix
//              X - vector rhs wrt continuous dofs
//
//  out:          
//              Y - result from multiplication wrt continous dofs
//
//  logic of the routine
//
//       X_DCG == CG_TO_DCG X         -> extend the continuous dofs 
//                                       to "discontinuous" dofs
//       Y_DCG == BlockHMat X_DCG     -> Y_DCG in "discontinous" format
//           Y == DCG_TO_CG Y_DCG     -> Y (squeezed Y_DCG) 
//
//   the rectangular connectivity matrices are given by the 
//   global arrays NR_CON and NODES_CON
//
// Global Variable used: 
// NDOFMACA, NDOFTOT (for testing)
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/GlobalVariables.h"

int zlusolC( complex double *y, complex double *x, iluptr lu );
void SolveTriDiag(long N,double complex *A,double complex *B,double complex *C,double complex *Y,double complex *F);

void BlockMHMatrixMultZ(double complex *X,double complex *Y,BlockHMatrix *BlockHMat)
{

  /* Local variables */
  int iprint=0;       // printing flag
  long N,M;           // dimension of MBH-Matrix
  long NB,MB;         // number of block MBH-matrices in row,column
  long NBlock,MBlock; // dimension of fixed Block
  long i,j,l,ib,jb;   // running indices  
  long IB,IE;         // ~
  long IPos,JPos;     // position of Block in HBH-Matrix
  long Error;
  double complex *XCopy,*YCopy; // copies of input taking into account
                                // the cluster permutation
  double complex *X_DCG,*Y_DCG; // MACA dof arrays

  HMatrix Block;                // local copy of a Block

  X_DCG=(double complex *)malloc(NDOFMACA*sizeof(double complex));
  assert(X_DCG != NULL);
  Y_DCG=(double complex *)malloc(NDOFMACA*sizeof(double complex));
  assert(Y_DCG != NULL);

  Expand_Vector(X, X_DCG);
  BlockHMatrixMult(*BlockHMat,X_DCG,Y_DCG);
  Squeeze_Vector(Y_DCG, Y);

  free(X_DCG);
  free(Y_DCG);

//  N=(*BlockHMat).NRow;
//  M=(*BlockHMat).NColumn;
//  
//  NB=(*BlockHMat).NBlockRow;
//  MB=(*BlockHMat).NBlockColumn;
//
//  Precond=(*BlockHMat).Precond; 
//
//  switch(Precond)
//  {
//    case None:
//      break;
//
//    case Jacobi:
//      for (i=0; i < N; i++)
//        Y[i]=Y[i]/(BlockHMat->PrecondDData)[i];
//      break;
//
//    case TriDiag:
//      SolveTriDiag(N-1,(BlockHMat->PrecondDData),(BlockHMat->PrecondDData)+N,(BlockHMat->PrecondDData)+2*N,Y,Y);       
//      break;
//
//    case ILU:
//      XCopy=(double complex *)malloc(N*sizeof(double complex));
//      assert(XCopy != NULL);
//
//      YCopy=(double complex *)malloc(N*sizeof(double complex));
//      assert(YCopy != NULL);
//
//      /* Permutation of the given vector with P_Row^T */
//      for (ib=0; ib < NB; ib++)
//	  {
//	    IPos  = (*BlockHMat).NDim[ib];
//	    Block = (*BlockHMat).HBlocks[ib];
//	    NBlock= Block.NRow;
//     
//        /* copy result take into account permutation of cluster points */
//        if(Block.NLocRow == NULL)
//        {
//          for (i=0; i < NBlock; i++)
//            YCopy[IPos+i]=Y[IPos+Block.PermuRow[i]];
//        }
//        else
//        {
//          if(iprint == 1) printf("BlockHMatrixMultZ\n");
//          l=0;
//          for (j=0; j < Block.NPermuRow; j++)
//          {
//            IB = Block.NLocRow[Block.PermuRow[j]];
//            IE = Block.NLocRow[Block.PermuRow[j]+1];
//            for(i=IB;i<IE;i++)
//            {
//              YCopy[IPos+l] = Y[IPos+i];
//              l++;
//              if(iprint == 1)
//              {
//                 printf("%4ld %4ld %4ld %14.8e + %14.8e*I\n",
//                        i,l,IPos,creal(YCopy[IPos+l]),cimag(YCopy[IPos+l]));
//                 fflush(stdout);
//              }
//            }
//          }
//          assert(l == NBlock);
//	    }
//	  }
//
//      /* Incomplete LU */
//      Error = zlusolC(YCopy,XCopy,BlockHMat->lu);
//      assert(Error == 0);
//
//      /* Permutation of the resulting vector with P_Column^T */
//      for (jb=0; jb < NB; jb++)
//	  {
//	    JPos  = (*BlockHMat).MDim[jb];
//	    Block = (*BlockHMat).HBlocks[jb*NB];
//	    MBlock= Block.NColumn;
//
//        /* copy result take into account permutation of cluster points */
//        if(Block.NLocColumn == NULL)
//        {
//          for (j=0; j < MBlock; j++)
//            Y[JPos+Block.PermuColumn[j]]=XCopy[JPos+j];
//        }
//        else
//        {
//          if(iprint == 1) printf("BlockHMatrixMultZ\n");
//          l=0;
//          for (j=0; j < Block.NPermuColumn; j++)
//          {
//            IB = Block.NLocColumn[Block.PermuColumn[j]];
//            IE = Block.NLocColumn[Block.PermuColumn[j]+1];
//            for(i=IB;i<IE;i++)
//            {
//              Y[JPos+i] = XCopy[JPos+l];
//              l++;
//              if(iprint == 1)
//              {
//                 printf("%4ld %4ld %4ld %14.8e + %14.8e*I\n",
//                        i,l,JPos,creal(XCopy[JPos+l]),cimag(XCopy[JPos+l]));
//                 fflush(stdout);
//              }
//            }
//          }
//          assert(l == MBlock);
//	    }
//	  }
//
//      /* Free memory */          
//      free(XCopy);
//      free(YCopy);
//      
//	  break;
//
//    case BasisElement:
//      XCopy=(double complex *)malloc(N*sizeof(double complex));
//      assert(XCopy != NULL);
//
//      cblas_zcopy(N,Y,1,XCopy,1);
//      SparseRowMult(N,(*BlockHMat).PrecondData,XCopy,Y);
//
//      /* Free memory */          
//      free(XCopy);
//      break;
//  }
}
