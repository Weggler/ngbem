/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                14.10.2011                                   *
\*****************************************************************************/
/*
  Last change 14.10.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void SparseRowMult(long N,SparseRow *SparseRows,double complex *x,double complex *y)
{
/*
  Local variables
*/
  long i,j;
  double complex C0=0.0+I*0.0;

  zset(N,C0,y,1);

  
/*   long *index; */
/*   FILE *Res; */
/*   char ResultName[4096]; */

/*   sprintf(ResultName,"Index%ld.dat",FekoProblem); */
/*   Res=fopen(ResultName,"w"); */

/*   index=(long *)malloc(N*sizeof(long)); */
/*   lset(N,0,index,1); */

/*   i=0; */
/*   for (j=0; j < SparseRows[i].NI; j++) */
/*     index[SparseRows[i].J[j]]=1; */

/*   fprintf(Res,"BEG_NDL_VARS\n"); */
/*   fprintf(Res,"index\n"); */
/*   for (i=0; i < N; i++) */
/*     fprintf(Res,"%ld\n",index[i]); */
/*   fprintf(Res,"END_NDL_VARS\n"); */

/*   fclose(Res); */

/*   exit(0); */

  for (i=0; i < N; i++)
    for (j=0; j < SparseRows[i].NI; j++)
      y[i]+=x[SparseRows[i].J[j]]*SparseRows[i].AIJ[j];
}
