/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                            Point-Line-Distance                              *
*                                                                             *
*                                10.12.2003                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "Includes.h"

double PLineDist(double X0[3], double X1[3], double X2[3])
{
    double X[3],Y[3],Z[3];

    cblas_dcopy(3,X2,1,X,1);
    cblas_daxpy(3,-1.0,X1,1,X,1);

    cblas_dcopy(3,X0,1,Y,1);
    cblas_daxpy(3,-1.0,X1,1,Y,1);

    dvecpr(X,Y,Z);

    return(cblas_dnrm2(3,Z,1)/cblas_dnrm2(3,X,1));

}
/*****************************************************************************\
*                                                                             *
*                                    The end                                  *
*                                                                             *
\*****************************************************************************/
