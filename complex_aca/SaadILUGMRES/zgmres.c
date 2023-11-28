#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "../CInclude/Includes.h"
#define  epsmac  1.0e-16

//#define dznrm2_ DZNRM2
//#define zdotc_ ZDOTC
//#define zaxpy_ ZAXPY
// ==begin== blas 1 prototypes
//double dznrm2_(int *, double complex *, int *);
//double complex zdotc_(int *, double complex *, int *, double complex *, int *);
//void zaxpy_(int *, double complex *, double complex *, int *, double complex *, int *);
// ====end== blas 1 prototypes

// ==begin== local prototypes
double sgn(double x, double y);
void zclartg(double complex f, double complex g, double *cs, double complex *sn,
             double complex *rot);
double abssq(double complex x);
void zgmres(int *pn, double complex *mat, double complex *rhs, double complex *sol, double *ptol, int *pim, int *pmaxits); 
// ====end== local prototypes

#define zgmres zgmres_
void zgmres(int *pn, double complex *mat, double complex *rhs, double complex *sol, double *ptol, int *pim, int *pmaxits){ 
/*----------------------------------------------------------------------
|                 *** GMRES calling zgemv ***                  
| on entry:
|========== 
|
| mat     = matrix 
|                                                                     
| rhs     = Complex vector of length n containing the right hand side.
|          
| sol     = Complex vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| pim     = Krylov subspace dimension 
| pmaxits = max number of iterations allowed. 
|
| on return:
|==========
| fgmresC  int =  0 --> successful return.
|          int =  1 --> convergence not achieved in itmax iterations.
|
| sol     == contains an approximate solution (upon successful return).
| pmaxits == has changed. It now contains the number of steps required
|            to converge -- 
+-----------------------------------------------------------------------
| work arrays: allocated internally and dynamically 
|=============
| vv    == work array of length [im+1][n] (used to store the Arnoldi
|          basis)
| hh    == work array of length [im][im+1] (Householder matrix)
| z     == work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
| armsol2, lusolD - preconditionning operation 
| BLAS1  routines.
|
+---------------------------------------------------------------------*/

  int n, maxits, im; 
  double tol;

  int i, i1, ii, j, k, k1, its, retval, tmp1 = 1;
  double complex **hh,  *s, *rs, t1, **vv, **z, negt, rot;
  double t, beta=1.0, eps1=0.0, *c; 
  double complex C1, C0;

  tol = *ptol;
  n = *pn;
  maxits = *pmaxits;
  im = *pim; 
  C1 = 1.0+0.0*I;
  C0 = 0.0+0.0*I;

//  printf("%d %d %d\n", *pn, *pim, *pmaxits);
//  for (j=0; j<1; j++)
//    printf("X: %d, % 25.14le + I %25.14le\n", j,creal(rhs[j]),cimag(rhs[j]));
//    for (i=0; i<n; i++)
//      printf("%d %d, % 25.14le + I %25.14le\n", i,j,creal(mat[i+j*n]),cimag(mat[i+j*n]));

   its = 0;
   vv = (double complex **) malloc((im+1)*sizeof(double complex *));
   z  = (double complex **) malloc(im*sizeof(double complex *));
   hh = (double complex **) malloc(im*sizeof(double complex *));
   for (i=0; i<im; i++) {
     vv[i] = NULL;
     z[i]  = NULL;
   }
   vv[im] = NULL;
   for (i=0; i<im; i++) {
     hh[i] = (double complex *) malloc((i+2)*sizeof(double complex));
   }
   c = (double *) malloc(im*sizeof(double ));
   s = (double complex *) malloc(im*sizeof(double complex));
   rs = (double complex *) malloc((im+1)*sizeof(double complex));
   
   
/*-------------------------------------------------------------
|   outer loop starts here
+------------------------------------------------------------*/
label20:
/*-------------------------------------------------------------
|   compute initial residual vector
+------------------------------------------------------------*/
   if( vv[0] == NULL )
       vv[0] = (double complex *)malloc( n*sizeof(double complex));
   cblas_zgemv(CblasColMajor,CblasNoTrans,n,n,
               &C1,mat,n,
               sol,1,&C0,vv[0],1);
   for (j=0; j<n; j++) {
     vv[0][j] = rhs[j] - vv[0][j];
   }
//   for (j=n-10; j<n; j++)
//     printf("%3d, %25.14le + I %25.14le\n", j, creal(vv[0][j]),cimag(vv[0][j]));
/*-------------------------------------------------------------
|  now vv[0] is the initial residual vector
+------------------------------------------------------------*/
//   beta = dznrm2_(&n, vv[0], &tmp1);
   beta = cblas_dznrm2(n, vv[0], tmp1);
//   printf("  Iter = %8d, Error = %15.8le\r",its,beta);
   if (beta == 0.0) goto label990;
   t = 1.0 / beta;
/*----------------------------------
|     normalize:  vv[0] = vv[0] / beta
+---------------------------------*/
   for (j=0; j<n; j++)
     vv[0][j] = vv[0][j]*t;
   if (its == 0) eps1 = tol*beta;
/*    ** initialize 1-st term  of rhs of hessenberg system   */
   rs[0] = beta + 0.0*I;
   i = -1;
 label4:
   i++;
   its++;
   i1 = i + 1;
/*------------------------------------------------------------
|    PRECONDITIONING    z_{j} = M^{-1} v_{j}
|                           w = A z_{j} = A M^{-1} v_{j}
+-----------------------------------------------------------*/
   if( z[i] == NULL )
     z[i] = (double complex *)malloc( n*sizeof(double complex));
   
   memcpy(z[i],vv[i],n*sizeof(double complex));

   if( vv[i1] == NULL) 
     vv[i1] = (double complex *)malloc( n*sizeof(double complex));
   cblas_zgemv(CblasColMajor,CblasNoTrans,n,n,
               &C1,mat,n,
               z[i],1,&C0,vv[i1],1);
//   for (j=0; j<10; j++)
//     printf("%3d, %25.14le + I %25.14le\n", j,creal(vv[i1][j]),cimag(vv[i1][j]));
//   for (j=n-10; j<n; j++)
//     printf("%3d, %25.14le + I %25.14le\n", j,creal(vv[i1][j]),cimag(vv[i1][j]));
//   printf("\n");
/*------------------------------------------------------------
|     modified gram - schmidt...
|     h_{i,j} = (w,v_{i})
|     w  = w - h_{i,j} v_{i}
+------------------------------------------------------------*/
   for (j=0; j<=i; j++) {
      //t1 = zdotc_(&n, vv[j], &tmp1, vv[i1], &tmp1);
      //zdotc_(&t1, &n, vv[j], &tmp1, vv[i1], &tmp1);
      cblas_zdotc_sub(n, vv[j], tmp1, vv[i1], tmp1, &t1);
      hh[i][j] = t1;
      negt = -t1;
      //zaxpy_(&n, &negt, vv[j], &tmp1, vv[i1], &tmp1);
      cblas_zaxpy(n, &negt, vv[j], tmp1, vv[i1], tmp1);
   }
/*----------------------------------
|     h_{j+1,j} = ||w||_{2}
+---------------------------------*/
//   t = dznrm2_(&n, vv[i1], &tmp1);
   t = cblas_dznrm2(n, vv[i1], tmp1);
   hh[i][i1] = t + 0.0*I;
   if (t == 0.0) goto label58;
   t = 1.0/t;
/*----------------------------------
|     v_{j+1} = w / h_{j+1,j}
+---------------------------------*/
   for (k=0; k<n; k++)
     vv[i1][k] = vv[i1][k]*t;
/*---------------------------------------------------
|     done with modified gram schimdt and arnoldi step
|     now  update factorization of hh
+--------------------------------------------------*/
 label58:
/*--------------------------------------------------------
|   perform previous transformations  on i-th column of h
+-------------------------------------------------------*/
   for (k=1; k<=i; k++) {
     k1 = k-1;
     t1 = hh[i][k1];
     hh[i][k1] = (c[k1])*t1 + s[k1]*hh[i][k];
     hh[i][k] = -conj(s[k1])*t1 + (c[k1])*hh[i][k];
    }
      
/*---------------------------------------------------
|     if gamma is zero then any small value will do...
|     will affect only residual estimate
+--------------------------------------------------*/
   //  if (gam == 0.0) gam = epsmac;
/*---------------------------------------------------
|     get  next plane rotation
+--------------------------------------------------*/

  zclartg(hh[i][i], hh[i][i1], &c[i], &s[i], &rot);
  // printf("c1 = %f \n", c[i]);
 
     rs[i1] = -conj(s[i])*rs[i];
     rs[i] =  (c[i])*rs[i];
/*----------------------------------------------------
|   determine residual norm and test for convergence
+---------------------------------------------------*/
//   hh[i][i] = conj(c[i])*hh[i][i] + s[i]*hh[i][i1];
	hh[i][i] = rot;
   beta = cabs(rs[i1]);
//   printf("  Iter = %8d, Error = %15.8le\r",its,beta);
   if ( (i < im-1) && (beta > eps1) && (its < maxits) )  goto label4;
/*---------------------------------------------------
|     now compute solution
|     first solve upper triangular system
+--------------------------------------------------*/
   rs[i] = rs[i]/hh[i][i];
   for (ii=1; ii<=i; ii++) {
     k=i-ii;
     k1 = k+1;
     t1=rs[k];
     for (j=k1; j<=i; j++){
       t1 = t1 - hh[j][k]*rs[j];
     	 rs[k] = t1/hh[k][k];
     	 }
   }
/*-----------------------------------------------------
|   form linear combination of v[i]'s to get solution   
+----------------------------------------------------*/
   for (j=0; j<=i; j++) {
     t1 = rs[j];
     for (k=0; k<n; k++)
       sol[k] += t1*z[j][k];
   }
//   for (k=0; k<n; k++)
//     printf("Y: %d, % 25.14le + I %25.14le\n", k, creal(sol[k]),cimag(sol[k]));
/*----------------------------------------------------
|     restart outer loop  when enecessary  
+---------------------------------------------------*/
   if (beta <= eps1) goto label990;
   if (its >= maxits) goto label991;
   goto label20;
label990:
   retval = 0;
   goto label888;
label991:
   retval = 1;
label888:
   for (i=0; i<=im; i++)
     if( vv[i] ) free(vv[i]);
   free(vv);
   for (i=0; i<im; i++){
     free(hh[i]);
     if( z[i] ) free(z[i]);
   }
   free(hh);
   free(z);
   free(c);
   free(s);
   free(rs);
   *pmaxits = its; 

   if(retval!=0)
   {
     printf("zgmres has not converged\n");
     fflush(stdout);
   }
//   else
//   {
//     printf("  Iter = %8d, Error = %15.8le\r",its,beta);
//     fflush(stdout);
//   }
}
/*-----------------end of gmres ------------------------------------
----------------------------------------------------------------------*/

