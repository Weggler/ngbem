//------------------------------------------------------------------
//
//  routine name  - GenMElement
//
//------------------------------------------------------------------
//
//  last revision - Feb 13
//
//  purpose       - interface to Fortran ...
//                  call elem_bem for mdle nodes
//                  ACA algorithm interface to Fortran ...
//
//  in:
//    IEl,JEl     - mdle node numbers
//    IDim,JDim   - dimension of the element matrix
//    
//  out:
//    Zloc        - element IEl-JEL single layer matrix 
//
//------------------------------------------------------------------

/* Includes for all functions */

# include "../CInclude/GlobalIncludes.h"

void maca_full_entry_(int *,int *,int *, double complex *);

void GenMElement(long IEl,long JEl,long IDim,long JDim, 
                 double complex *Zloc)
{
/* 
  Local Variables
*/
  int i, j;
  int maxdof;
  double complex *zsloc;
  if(IDim < JDim) 
    maxdof = (int) JDim;
  else
    maxdof = (int) IDim;
  zsloc = (double complex*) malloc(maxdof*maxdof*sizeof(double complex));

/* 
  1) C-numbering to Fortran-numbering
  2) generate element matrices 
  3) copy result column-wise
*/
  i = 1+(int) IEl;
  j = 1+(int) JEl;

  maca_full_entry_(&i,&j,&maxdof, zsloc);

  for(j=0;j<JDim;j++) 
  {
    for(i=0;i<IDim;i++)
    {
      Zloc[i+j*JDim] = zsloc[i+j*JDim]; 
//      printf("%ld %ld %lf %lf\n", IEl,JEl, __real__ zsloc[i+j*JDim], __imag__ zsloc[i+j*JDim]);
    }
  }
}
