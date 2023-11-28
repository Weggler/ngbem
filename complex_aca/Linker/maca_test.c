//------------------------------------------------------------------
//
//  routine name  -  maca_tset
//
//------------------------------------------------------------------
//
//  last revision -  May 12
//
//  purpose       -  compute || A_c - P^t A_d P ||_F, where
//                    A_c == system matrix contin. Galerkin
//                    A_d == system matrix discontin. Galerkin
//                    P,P^t == connectivity matrices
//                   
//
//  in:
//           ncol - number of column to compare 
//         zbglob - ncol column of A_c
//        ndoftot - number of total degrees of freedom
//
//  out:
//         zbglob - difference
//         norm   - contribution of this column to Frobenius norm
//
//------------------------------------------------------------------

/* Includes global Variables (System Matrix!) */
# include "../CInclude/GlobalVariables.h"

#define maca_test maca_test_
void maca_test(int *ncol, double complex *zbglob, int *ndoftot, double *norm)
{
  int iprint = 0;
  int nn;
  long n,i,j;
  double complex *zsglob, *test;
  double cnorm = 0.0;
 
  j  = *ncol;
  n  = (long) *ndoftot;
  nn = *ndoftot;
 
  zsglob = (double complex*) malloc(n*sizeof(double complex));
  test = (double complex*) malloc(n*sizeof(double complex));

  zset(nn,0.0+I*0.0,zsglob,1);
  zset(nn,0.0+I*0.0,test,1);
  zsglob[j-1] = 1.0 + I*0.0;

  BlockMHMatrixMultZ(zsglob,test,BHMAT);

  for(i=0;i<n;i++) 
  {
    cnorm += (test[i]-zbglob[i])*conj(test[i]-zbglob[i]);
  }
//  printf("maca_test: %d %25.14le\n", j, cnorm);
  *norm = cnorm;

  free(zsglob);
  free(test);
}
