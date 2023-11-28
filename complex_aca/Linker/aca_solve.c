//------------------------------------------------------------------
//
//  routine name  -  aca_solve
//
//------------------------------------------------------------------
//
//  last revision -  Aug 12
//
//  purpose       -  solve the linear system by
//                   BlockHLUSolveMM(...) (direct    solver)
//                or zfgmres1(...)        (iterative solver)
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

#define aca_solve aca_solve_
void aca_solve(double complex *zbglob, int *ndoftot)
{
  int nn;
  int GMResMaxIter, Result;
  long n,i,j;
  double ACAEps, GMResEps;
  double complex *zsglob;

  n  = (long) *ndoftot;
  nn = *ndoftot;
  ACAEps = PARAMACA.ACAEps;

  zsglob = (double complex*) malloc(n*sizeof(double complex));

///* Direct solver */
//  printf("Start Block LU Solver ...\n"); 
//  fflush(stdout);
//
//  BlockHLUSolveMM(1,BHMAT,zsglob,zbglob,ACAEps);
//  for(i=0;i<n;i++)
//    zbglob[i] = zsglob[i];
//
//  printf("  Done ...\n"); 
//  fflush(stdout);

/* GMRes solver without preconditioner*/
  printf("Start GMRes Solver ...\n"); 
  fflush(stdout);

  GMResEps    =1.0e-8;
  GMResMaxIter=20000;

  zset(nn,0.0+I*0.0,zsglob,1);

  Result = zfgmres1(nn,(void (*)(double complex *,double complex *,void *))BlockHMatrixMultZ,BHMAT,zbglob,zsglob,GMResEps,GMResMaxIter,&GMResMaxIter,GMResPrint);

  for(i=0;i<n;i++)
    zbglob[i] = zsglob[i];

  printf("  \n");
  printf("  Done with %d iterations\n\n",GMResMaxIter);
  fflush(stdout); 

///* Preconditioned GMRes solver */
//  printf("Start GMRes Solver with Preconditioning ...\n"); 
//  (*BHMAT).Precond = BasisElement;
//  fflush(stdout);
//
//  BlockPrePrecond(BHMAT,GenElement, zbglob);
//
//  GMResEps    =1.0e-8;
//  GMResMaxIter=20000;
//
//  zset(nn,0.0+I*0.0,zsglob,1);
//
//  Result = zfgmres1(nn,(void (*)(double complex *,double complex *,void *))BlockHMatrixMultZ,BHMAT,zbglob,zsglob,GMResEps,GMResMaxIter,&GMResMaxIter,GMResPrint);
//
//  for(i=0;i<n;i++)
//    zbglob[i] = zsglob[i];
//
//  printf("  \n");
//  printf("  Done with %d iterations\n\n",GMResMaxIter);
//  fflush(stdout); 

  /* free local variables */
  free(zsglob);
}
