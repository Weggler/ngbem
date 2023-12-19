/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                07.06.2010                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/


double complex Element11(long IEl,long JEl)
{
/*
  Feko variables
*/
  integer num_m,num_n;
  Fekodoublecomplex ctemp;

  num_m=(integer) IEl+1;
  num_n=(integer) JEl+1;

  feko_get_amn__ (&num_m, &num_n, &ctemp);

  return ctemp.r+I*ctemp.i;
}
