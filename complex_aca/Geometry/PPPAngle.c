/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                         Point-Point-Point-Angle                             *
*                                                                             *
*                                31.12.2003                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "Includes.h"

double PPPAngle(double X[3], double Y[3], double Z[3])
{
    double X1[3],X2[3];

    cblas_dcopy(3,X,1,X1,1);
    cblas_daxpy(3,-1.0,Y,1,X1,1);

    cblas_dcopy(3,Z,1,X2,1);
    cblas_daxpy(3,-1.0,Y,1,X2,1);

    return(acos(cblas_ddot(3,X1,1,X2,1)/(cblas_dnrm2(3,X1,1)*cblas_dnrm2(3,X2,1))));

}
/*****************************************************************************\
*                                                                             *
*                                    The end                                  *
*                                                                             *
\*****************************************************************************/
