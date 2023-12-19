//------------------------------------------------------------------
//                                                                  
//  routine name  -  Expand_Vector, Squeeze_Vector
//                                                                  
//------------------------------------------------------------------
//
//  last revision -  May 12
//
//  purpose       -  expands (squeezes) the given vector by global 
//                   connectivity data. 
//
//  in/out:
//              X - vector rhs wrt continuous Galerkin dofs 
//          X_DCG - vector rhs wrt discontinuous Galerkin dofs == MACA dofs 
//
//
//  Global Variables: 
//        NDOFTOT - nr of cG dof's
//       NDOFMACA - nr of dcG dof's
//         NR_CON - array of length NDOFTOT nr of dcG dofs per cG dofs  
//      NODES_CON - contains the MACA dcG dof numbers that belong to
//                  the cG dofs
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/GlobalVariables.h"

void Expand_Vector(double complex *X, double complex *X_DCG)
{
  /* local variables */
  int i=0;
  int k,l,kE;

  for(k=0;k<NDOFTOT;k++)
  {
    for(l=0;l<NR_CON[k];l++)
    {
      kE = NODES_CON[i];
      X_DCG[kE] = X[k];
      i++;
    }
  }
}

void Squeeze_Vector(double complex *X_DCG, double complex *X)
{
  /* local variables */
  int i=0;
  int l,k,kE;

  for(k=0;k<NDOFTOT;k++)
  {
    X[k] = 0.0;
    for(l=0;l<NR_CON[k];l++)
    {
      kE = NODES_CON[i];
      X[k] = X[k]+X_DCG[kE];
//      if(k==0)
//        printf("Squeeze Vector %25.14le + I %25.14le\n",
//               creal(X_DCG[kE]), cimag(X_DCG[kE]));
      i++;
    }
  }
}
