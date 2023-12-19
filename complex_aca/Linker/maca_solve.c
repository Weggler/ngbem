//------------------------------------------------------------------
//
//  routine name  - maca_solve
//
//------------------------------------------------------------------
//
//  last revision - Feb 13
//
//  purpose       - solve the linear system by iterative
//                  zfgmres1(...)  solver
//
//  in:
//         zbglob - rhs
//        ndoftot - number of total degrees of freedom
//
//  out:
//         zbglob - solution
//
//------------------------------------------------------------------

/* Includes global Variables (System Matrix!) */
# include "../CInclude/GlobalVariables.h"

#define maca_solve maca_solve_
void maca_solve(double complex *zbglob, int *ndoftot)
{
  int iprint = 0;
  int nn;
  int GMResMaxIter, Result;
  long n,i,j;
  double ACAEps, GMResEps;
  double complex *zsglob;
 
  n  = (long) *ndoftot;
  nn = *ndoftot;
  ACAEps = PARAMACA.ACAEps;
 
  zsglob = (double complex*) malloc(n*sizeof(double complex));

/* GMRes solver without preconditioner*/
  printf("Start GMRes Solver ...\n"); 
  fflush(stdout);

  GMResEps    =1.0e-6;
  GMResMaxIter=20000;
  printf("maca_init: the GMRes parameters are ...\n");
  printf("           GMResEps     = %le\n", GMResEps);
  printf("           GMResMaxIter = %d\n",  GMResMaxIter);
 
  zset(nn,0.0+I*0.0,zsglob,1);
  Result = zfgmres1(nn,
           (void (*)(double complex *,double complex *,void *))BlockMHMatrixMultZ,
                    BHMAT,zbglob,zsglob,
                    GMResEps,GMResMaxIter,&GMResMaxIter,GMResPrint);
  assert(Result == 0);

  for(i=0;i<n;i++)
    zbglob[i] = zsglob[i];

  if(iprint==1)
  {
    for(i=0;i<n;i++)
      printf("RESULT: %4ld % 6.5lf+ I % 6.5lf\n",
             i,creal(zbglob[i]),cimag(zbglob[i]));
  }

  printf("  \n");
  printf("  Done with %d iterations\n\n",GMResMaxIter);
  fflush(stdout); 

  /* free local variables */
  free(zsglob);
}
