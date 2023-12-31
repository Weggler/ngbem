/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                          Laplace equation in 3D                             *
*                                                                             *
*                                28.12.2003                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "Includes.h"

void dvecpr(double X[3], double Y[3], double Z[3])
{
    Z[0]=X[1]*Y[2]-X[2]*Y[1];
    Z[1]=X[2]*Y[0]-X[0]*Y[2];
    Z[2]=X[0]*Y[1]-X[1]*Y[0];

}
/*****************************************************************************\
*                                                                             *
*                                    The end                                  *
*                                                                             *
\*****************************************************************************/
