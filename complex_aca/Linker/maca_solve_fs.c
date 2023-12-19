//------------------------------------------------------------------
//
//  routine name  -  maca_solve_fs
//
//------------------------------------------------------------------
//
//  last revision -  Aug 12
//
//  purpose       -  solve the linear system by iterative
//                   zfgmres1(...)  solver
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

#define maca_solve_fs maca_solve_fs_
void maca_solve_fs(double complex *zbglob, int *ndoftot)
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

//  for (i=0; i < 20; i++)
//    printf("i,x[i]= %ld %25.15e %25.15e\n",i,creal(zbglob[i]),cimag(zbglob[i]));
//  printf("\n");
//  for (i=nn-20; i<nn; i++) 
//    printf("i,x[i]= %ld %25.15e %25.15e\n",i,creal(zbglob[i]),cimag(zbglob[i]));

//  /* compute the RHS */
//  for(i=0;i<n;i++)
//  {
//    zsglob[i] = zbglob[i];
//    zbglob[i] = 0.0 + I*0.0;
//  }
//  BlockHMatrixMult(*BHMAT,zsglob,zbglob);
//  if(iprint==1)
//  {
//    for(i=0;i<n;i++)
//      printf("RHS: %4ld % 6.5lf+ I % 6.5lf\n",
//             i,creal(zbglob[i]),cimag(zbglob[i]));
//  }

///* GMRes solver without preconditioner*/
//  printf("Start GMRes Solver ...\n"); 
//  fflush(stdout);
//
//  GMResEps    =1.0e-16;
//  GMResMaxIter=20000;
//
//  zset(nn,0.0+I*0.0,zsglob,1);
//  Result = zfgmres1(nn,
//           (void (*)(double complex *,double complex *,void *))BlockMHMatrixMultZ_fs,
//                    BHMAT,zbglob,zsglob,
//                    GMResEps,GMResMaxIter,&GMResMaxIter,GMResPrint);
//  assert(Result == 0);
//
//  printf("  \n");
//  printf("  Done with %d iterations\n\n",GMResMaxIter);
//  fflush(stdout); 
//
//  /* print the first 20 and the last 20 entries in result vector... 
//     they decrease and thus, it's good to check on how many relevant
//     positions after the decimal point are available */
//  for (i=0; i < 20; i++)
//    printf("i,x[i]= %ld %25.15e %25.15e\n",i,creal(zsglob[i]),cimag(zsglob[i]));
//  printf("\n");
//  for (i=nn-20; i<nn; i++) 
//    printf("i,x[i]= %ld %25.15e %25.15e\n",i,creal(zsglob[i]),cimag(zsglob[i]));

  long ipos,jpos,II,JJ,N; 
  double complex fronorm,cerror; 
  double complex A[9], aij; 
 
  cerror =0.0+0.0*I; 
  fronorm=0.0+0.0*I; 

  N = nn/3;
  /* compute the real error of the approximation */
  for (JJ=0; JJ<N; JJ++) 
  { 
    for (j=0; j<3; j++) 
    { 
      jpos = 3*JJ+j;
      zset(nn,0.0+0.0*I,zbglob,1); 
      zbglob[jpos] = 1.0+0.0*I; 
      BlockHMatrixMult(*BHMAT,zbglob,zsglob);
      for (II=0; II<N; II++) 
      { 
        GenMElement(II,JJ,3,3,A); 
        for (i=0; i<3; i++) 
        { 
          ipos = 3*II+i;
          aij = A[i+3*j];
          fronorm += aij*conj(aij); 
          cerror  += (zsglob[ipos]-aij)*conj(zsglob[ipos]-aij); 
        } 
      } 
    } 
  } 
  printf("Fro-Norm = %25.15e Rel-Error = %25.15e\n",
         sqrt(creal(fronorm)),sqrt(creal(cerror/fronorm))); 

//  /* check single entries */
//  while(1 == 1) 
//  { 
// 	printf("i="); 
// 	scanf("%ld",&ipos); 
//      if (ipos == -1) break; 
// 	printf("j="); 
// 	scanf("%ld",&jpos); 
//      
//    zset(nn,0.0+0.0*I,zbglob,1); 
//    zbglob[jpos]=1.0+0.0*I; 
//
//    BlockHMatrixMult(*BHMAT,zbglob,zsglob);
//    II = ipos / 3; i = ipos % 3;
//    JJ = jpos / 3; j = jpos % 3;
//    GenMElement(II,JJ,3,3,A); 
//    aij = A[i+3*j];
//    printf("App = %14.8e + %14.8e*I\nExa = %14.8e + %14.8e*I\n",
//	       creal(zsglob[ipos]),cimag(zsglob[ipos]),
//           creal(aij),cimag(aij)); 
//  } 

  for(i=0;i<n;i++)
    zbglob[i] = zsglob[i];

  /* free local variables */
  free(zsglob);
}
