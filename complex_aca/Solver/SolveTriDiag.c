/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                09.05.2010                                   *
\*****************************************************************************/

# include "../CInclude/Includes.h"

void SolveTriDiag(long N,double complex *A,double complex *B,double complex *C,double complex *Y,double complex *F)
{
  double complex *al,*be;
  int i;
   
  al = (double complex *)malloc((N+1)*sizeof(double complex));
  assert(al != NULL);
  be = (double complex *)malloc((N+1)*sizeof(double complex));
  assert(be != NULL);
   
  al[0]=B[0]/C[0]; 
  be[0]=F[0]/C[0];

  for (i = 1; i < N; i++)
    {
      al[i]=B[i]/(C[i]-A[i]*al[i-1]);
      if (creal(al[i])*creal(al[i])+cimag(al[i])*cimag(al[i]) > 1.0) 
	al[i]=1.0+0.0*I;
      be[i]=(F[i]+A[i]*be[i-1])/(C[i]-A[i]*al[i-1]);
/*       printf("ai= %le %le\n",creal(A[i]),cimag(A[i])); */
/*       printf("bi= %le %le\n",creal(B[i]),cimag(B[i])); */
/*       printf("ci= %le %le\n",creal(C[i]),cimag(C[i])); */
/*   printf("al= %le %le\n",creal(al[i]),cimag(al[i])); */
/*   printf("ci= %le %le\n",creal(C[i]-A[i]*al[i-1]),cimag(C[i]-A[i]*al[i-1])); */
/*   getchar(); */
    }

  Y[N]=(F[N]+A[N]*be[N-1])/(C[N]-A[N]*al[N-1]);

  for (i = N-1; i >= 0; i--)
    {
      Y[i]=al[i]*Y[i+1]+be[i];
    }

  free(al); free(be); 
}
