//------------------------------------------------------------------
//                                                                  
//  routine name  -  BlockHMatrixMult
//                                                                  
//------------------------------------------------------------------
//
//  last revision -  Feb 13
//
//  purpose       -  matrix vector multiplication of a
//                   Block MH-matrix with a vector
//
//  in:
//           HMat - BMH-Matrix
//              X - vector
//
//  out:          
//              Y == HMat X
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

void BlockHMatrixMult(BlockHMatrix BlockHMat,double complex *X,double complex *Y)
{
  /* Local variables */
  int iprint=0;
  long N,M;
  long NB,MB;
  long ib,jb;
  long Block,XPos,YPos,MY;
  double complex *YCopy;
  const double complex One=1.0+I*0.0;

  /* Initialisation */
  N=BlockHMat.NRow;
  M=BlockHMat.NColumn;
  
  NB=BlockHMat.NBlockRow;
  MB=BlockHMat.NBlockColumn;

  zset(N,0.0+I*0.0,Y,1);

  Block=0;

  
  for (jb=0; jb<MB; jb++)
  {
    XPos=BlockHMat.NDim[jb];

    if(iprint == 1)
    {
      printf("%s X: %4ld % 6.5lf+ I % 6.5lf\n",
             __func__, jb,creal(X[XPos+jb]),cimag(X[XPos+jb]));
    }
  
    for (ib=0; ib < NB; ib++)
    {
       YPos=BlockHMat.NDim[ib];
  
       MY=(BlockHMat.HBlocks[Block]).NRow;
  
       YCopy=(double complex *)malloc(MY*sizeof(double complex));
       assert(YCopy != NULL);
  
       HMatrixMult(BlockHMat.HBlocks[Block],X+XPos,YCopy);
       cblas_zaxpy(MY,&One,YCopy,1,Y+YPos,1);
  
       free(YCopy);
  
       Block++;    
    }  
  }
  if(iprint == 1)
  {
     for(ib=0; ib<MY; ib++)
       printf("%s Y: %4ld % 6.5lf+ I % 6.5lf\n",
              __func__,ib,creal(Y[YPos+ib]),cimag(Y[YPos+ib]));
  }
  return;
}
