/**********************************************************************\
*                                                                      *
*                          S. Rjasanow C-Software                      *
*                                                                      *
*                       Adaptive Cross Approximation                   *
*                                                                      *
*                                Mar 2012                              *
\**********************************************************************/

/* Includes for all functions */

# include "../CInclude/GlobalIncludes.h"

void zaca_full_entry_(int *, int *, double complex *);

double complex GenElement(long IEl,long JEl)
{
  int i;
  int j;
  double complex zvalue;

// C-numbering to Fortran-numbering
  i = 1+(int) IEl;
  j = 1+(int) JEl;
  zaca_full_entry_(&i, &j, &zvalue);
//  if(IEl == 0)
//  {
//    printf("%ld %ld\n", IEl, JEl);
//    printf("%25.16lf %25.16lf\n", __real__ zvalue, __imag__ zvalue);
//  }
//  if(i==j)
//    zvalue = 1.;
//  else 
//    zvalue = 0.;

  return zvalue;
}
