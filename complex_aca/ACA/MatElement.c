/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                08.03.2009                                   *
\*****************************************************************************/
/*
  Last change 08.03.2009
*/

/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

double complex MatElement(long NBlock,long MBlock,long Rank,double complex *U,double complex *VT,long i,long j)
{
/*
  Local variables
*/
  long k;
  double complex Sum;

  Sum=0.0+0.0*I;

  if (Rank > 0)
    for (k=0; k < Rank; k++)
      Sum+=U[k*NBlock+i]*conj(VT[k*MBlock+j]);

  return Sum;
}

//double complex MatMElement(long NBlock,long MBlock,long Rank,double complex *U,double complex *VT,long II,long JJ,long ii,long jj)
//{
///*
//  Local variables
//*/
//  long k;
//  double complex Sum;
//
//  Sum=0.0+0.0*I;
//
//  if (Rank > 0)
//    for (k=0; k < Rank; k++)
//    {
//      Sum+=U[k*NBlock+II+ii]*conj(VT[(k+II)*MBlock+JJ+jj]);
//      printf("%ld:\t % lf +I % lf\n", (k+II)*MBlock+JJ+jj, 
//      creal(VT[(k+II)*MBlock+JJ+jj]), 
//      cimag(VT[(k+II)*MBlock+JJ+jj]));
//    }
//
//  return Sum;
//}
