// the solver
int zfgmres1(int nn, void (*matvec)(double complex *,double complex *,void*), void* param,  double complex *rhs, 
	    double complex *sol, double tol, int im, int *itmax,  void (*prn_fnc)(int,double)) ;
/*----------------------------------------------------------------------
|                 *** Preconditioned FGMRES ***                  
+-----------------------------------------------------------------------
| This is a simple Complex version of the ARMS preconditioned FGMRES algorithm. 
+-----------------------------------------------------------------------
| on entry:
|========== 
|
|(matvec) = matrix vector multiplication routine
| matvec(void* param, double complex *in, double complex *out);
|
| param   = eventual 
|                                                                      
| rhs     = Complex vector of length n containing the right hand side.
|           
| sol     = Complex vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| im      = Krylov subspace dimension 
| itmax   = max number of iterations allowed. 
|
| on return:
|==========
| fgmresC  int =  0 --> successful return.
|          int =  1 --> convergence not achieved in itmax iterations.
|
| sol   == contains an approximate solution (upon successful return).
| itmax == has changed. It now contains the number of steps required
|          to converge -- 
+-----------------------------------------------------------------------
| work arrays: allocated internally and dynamically 
|=============
| vv    == work array of length [im+1][n] (used to store the Arnoldi
|          basis)
| hh    == work array of length [im][im+1] (Householder matrix)
| z     == work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
| armsol2, lusolD - preconditionning operation 
| BLAS1  routines.
|
+---------------------------------------------------------------------*/
