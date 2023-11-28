#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
//#include "mkl_spblas.h"
#include "mkl_service.h"
//#include "iccsolvers.h"


int solve_system_gmres_mv(long nEq, 
		void (*mat_vec)(double*,double*,void*), void* data,
		void (*pre_con)(double*,double*,void*), void* pdata,
		double* rhs, double* solution,
		double eps, long maxIter,
		void (*prn_fnc)(int,double)) {

	/*---------------------------------------------------------------------------
	/ Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	/---------------------------------------------------------------------------*/
	MKL_INT ipar[128];
	double dpar[128];
	long tmp_size= ((2*maxIter+1)*nEq + maxIter*(maxIter+9)/2 + 1);
	double *sol, *tmp, *residual;
	/*---------------------------------------------------------------------------
	/ Some additional variables to use with the RCI (P)FGMRES solver
	/---------------------------------------------------------------------------*/
	MKL_INT itercount;
	MKL_INT RCI_request, i, ivar;
	double dvar;
	//printf("--------------------------------------------------------\n");
	//printf("The FULLY ADVANCED example of usage of RCI FGMRES solver\n");
	//printf("   to solve a non-symmetric indefinite non-degenerate\n");
	//printf("          algebraic system of linear equations\n");
	//printf("--------------------------------------------------------\n\n");
	//alloc the tmp
	tmp      = (double*) malloc(tmp_size*sizeof(double));
	sol      = (double*) malloc(     nEq*sizeof(double));
	residual = (double*) malloc(     nEq*sizeof(double));
	/*---------------------------------------------------------------------------
	/ Initialize variables and the right hand side through matrix-vector product
	/---------------------------------------------------------------------------*/
	ivar=nEq;
	/*---------------------------------------------------------------------------
	/ Initialize the solver
	/---------------------------------------------------------------------------*/
	dfgmres_init(&ivar, solution, rhs, &RCI_request, ipar, dpar, tmp);
	if (RCI_request!=0) goto FAILED;
	/*---------------------------------------------------------------------------
	/  Set the desired parameters:
	/---------------------------------------------------------------------------*/
	ipar[14]=maxIter; // do the restart after maxIter iterations
	ipar[4]=maxIter; // set max iteration count to maxIter
	ipar[7]=1;       // do the stopping based on max iter count
	if(pre_con!=NULL)
		ipar[10]=0;      // use the non-preconditioned version of the gmres
	else
		ipar[10]=1;      // use the preconditioned version of the gmres
	dpar[0]=eps;    // set the target relative tolerance to eps
	/*---------------------------------------------------------------------------
	/  Check the correctness and consistency of the newly set parameters
	/ ---------------------------------------------------------------------------*/
	dfgmres_check(&ivar, solution, rhs, &RCI_request, ipar, dpar, tmp);
	if (RCI_request!=0) goto FAILED;
	/*---------------------------------------------------------------------------
	/  Compute the solution by RCI (P)FGMRES solver with preconditioning
	/  Reverse Communication starts here
	/ ---------------------------------------------------------------------------*/
ONE:  dfgmres(&ivar, solution, rhs, &RCI_request, ipar, dpar, tmp);
	/*---------------------------------------------------------------------------
	/  If RCI_request=0, then the solution was found with the required precision
	/ ---------------------------------------------------------------------------*/
	if (RCI_request==0) goto COMPLETE;
	/*---------------------------------------------------------------------------
	/  If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/  and put the result in vector tmp[ipar[22]-1]
	/ ---------------------------------------------------------------------------
	/  NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
	/  therefore, in C code it is required to subtract 1 from them to get C style
	/  addresses
	/ ---------------------------------------------------------------------------*/
	if (RCI_request==1)
	{
		//mkl_dcsrgemv(&cvar, &ivar, A, ia, ja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
		mat_vec(&tmp[ipar[21]-1],&tmp[ipar[22]-1],data);
		goto ONE;
	}
	/*---------------------------------------------------------------------------
	/  If RCI_request=2, then do the user-defined stopping test
	/  The residual stopping test for the computed solution is performed here
	/ ---------------------------------------------------------------------------*/
	if (RCI_request==2)
	{
		if(prn_fnc!=NULL) prn_fnc(ipar[3],dpar[4]);
		if (dpar[4]<eps) goto COMPLETE;
		else goto ONE;
		/* Request to the dfgmres_get routine to put the solution into sol[N] via ipar[12]
		--------------------------------------------------------------------------------
		WARNING: beware that the call to dfgmres_get routine with ipar[12]=0 at this
		stage may destroy the convergence of the FGMRES method, therefore, only
		advanced users should exploit this option with care */
		ipar[12]=1;
		/* Get the current FGMRES solution in the vector sol */
		dfgmres_get(&ivar, solution, sol, &RCI_request, ipar, dpar, tmp, &itercount);
		/* Compute the current true residual via MKL (Sparse) BLAS routines */
		//mkl_dcsrgemv(&cvar, &ivar, A, ia, ja, sol, residual);
		mat_vec(sol,residual,data);
		dvar= - 1.0E0;
		i=1;
		daxpy(&ivar, &dvar, rhs, &i, residual, &i);
		dvar=dnrm2(&ivar,residual,&i);
		if(prn_fnc!=NULL) prn_fnc(itercount,dvar);
		if (dvar<eps) goto COMPLETE;
		else goto ONE;
	}
	/*---------------------------------------------------------------------------
	/  If RCI_request=3, then apply the preconditioner on the vector
	/  tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]
	/ ---------------------------------------------------------------------------
	/  NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
	/  therefore, in C code it is required to subtract 1 from them to get C style
	/  addresses
	/ ---------------------------------------------------------------------------*/
	if (RCI_request==3)
	{
		if(pre_con!=NULL) {
			pre_con(&tmp[ipar[21]-1],&tmp[ipar[22]-1],pdata);
		} else {
			i=1;
			dcopy(&ivar, &tmp[ipar[21]-1], &i, &tmp[ipar[22]-1], &i);
		}
		goto ONE;
	}
	/*---------------------------------------------------------------------------
	/  If RCI_request=4, then check if the norm of the next generated vector is
	/  not zero up to rounding and computational errors. The norm is contained
	/  in dpar[6] parameter
	/ ---------------------------------------------------------------------------*/
	if (RCI_request==4)
	{
		if (dpar[6]<1.0E-12) goto COMPLETE;
		else goto ONE;
	}
	/*---------------------------------------------------------------------------
	/  If RCI_request=anything else, then dfgmres subroutine failed
	/  to compute the solution vector: computed_solution[N]
	/ ---------------------------------------------------------------------------*/
	else
	{
		goto FAILED;
	}
	/*---------------------------------------------------------------------------
	/  Reverse Communication ends here
	/  Get the current iteration number and the FGMRES solution (DO NOT FORGET to
	/  call dfgmres_get routine as computed_solution is still containing
	/  the initial guess!). Request to dfgmres_get to put the solution
	/  into vector computed_solution[N] via ipar[12]
	/ ---------------------------------------------------------------------------*/
COMPLETE:  
	ipar[12]=0;
	dfgmres_get(&ivar, solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
	/*-------------------------------------------------------------------------*/
	/* Release internal MKL memory that might be used for computations         */
	/* NOTE: It is important to call the routine below to avoid memory leaks   */
	/* unless you disable MKL Memory Manager                                   */
	/*-------------------------------------------------------------------------*/
	free(sol);
	free(tmp);
	free(residual);
        mkl_free_buffers();
	return itercount;

	/*-------------------------------------------------------------------------*/
	/* Release internal MKL memory that might be used for computations         */
	/* NOTE: It is important to call the routine below to avoid memory leaks   */
	/* unless you disable MKL Memory Manager                                   */
	/*-------------------------------------------------------------------------*/
FAILED: fprintf(stderr,"solve_system_gmres_mv has FAILED as the solver has returned the ERROR ");
	fprintf(stderr,"code %d", RCI_request);
	free(sol);
	free(tmp);
	free(residual);
        mkl_free_buffers();
        return  - 1;
}

/*---------------------------------------------------------------------------*/
/*  Function for solving symmetric positive definite system of equations.    */
/*  Simplest case: no preconditioning and no the user-defined stopping tests.*/
/*---------------------------------------------------------------------------*/
int solve_sym_system_cg_mv(long nEq, void (*mv)(double*,double*,void*), void* data, double* rhs, double* solution, double eps, long maxIter) {
	/*---------------------------------------------------------------------------*/
	/* Define arrays for the upper triangle of the coefficient matrix and rhs vector*/
	/* Compressed sparse row storage is used for sparse representation           */
	/*---------------------------------------------------------------------------*/
	MKL_INT n=(MKL_INT)nEq, rci_request, itercount;//, expected_itercount=8, i;
	/* Fill all arrays containing matrix data. */

	/*---------------------------------------------------------------------------*/
	/* Allocate storage for the solver ?par and temporary storage tmp            */
	/*---------------------------------------------------------------------------*/
	//MKL_INT length=128;
	/*---------------------------------------------------------------------------*/
	/* Some additional variables to use with the RCI (P)CG solver                */
	/*---------------------------------------------------------------------------*/
	MKL_INT ipar[128];
	double  dpar[128],*tmp;
	//double eone=-1.E0;
	char tr='u';
	//MKL_INT ione=1;
	
	tmp = (double*) malloc(4*n*sizeof(double));
	/*---------------------------------------------------------------------------*/
	/* Initialize the solver                                                     */
	/*---------------------------------------------------------------------------*/
	dcg_init(&n,solution,rhs,&rci_request,ipar,dpar,tmp);
	//fprintf(stderr,"code %d\n", rci_request);
	if (rci_request!=0) goto failure;
	/*---------------------------------------------------------------------------*/
	/* Set the desired parameters:                                               */
	/* LOGICAL parameters:                                                       */
	/* do residual stopping test                                                 */
	/* do not request for the user defined stopping test                         */
	/* DOUBLE parameters                                                         */
	/* set the relative tolerance to eps    instead of default value 1.0D-6      */
	/*---------------------------------------------------------------------------*/
	ipar[4]=maxIter; //max iter = maxIter
	ipar[8]=1;
	ipar[9]=0;
	dpar[0]=eps;
	/*---------------------------------------------------------------------------*/
	/* Check the correctness and consistency of the newly set parameters         */
	/*---------------------------------------------------------------------------*/
	dcg_check(&n,solution,rhs,&rci_request,ipar,dpar,tmp);
	//fprintf(stderr,"code %d\n", rci_request);
	if (rci_request!=0) goto failure;
	/*---------------------------------------------------------------------------*/
	/* Compute the solution by RCI (P)CG solver without preconditioning          */
	/* Reverse Communications starts here                                        */
	/*---------------------------------------------------------------------------*/
rci: dcg(&n,solution,rhs,&rci_request,ipar,dpar,tmp);
	//fprintf(stderr,"code %d\n", rci_request);
	/*---------------------------------------------------------------------------*/
	/* If rci_request=0, then the solution was found with the required precision */
	/*---------------------------------------------------------------------------*/
	if (rci_request==0) goto getsln;
	/*---------------------------------------------------------------------------*/
	/* If rci_request=1, then compute the vector A*tmp[0]                        */
	/* and put the result in vector tmp[n]                                       */
	/*---------------------------------------------------------------------------*/
	if (rci_request==1)
	{
		mv(tmp,&(tmp[n]),data);
		//mkl_dcsrsymv(&tr, &n, a, ia, ja, tmp, &tmp[n]);
		goto rci;
	}
	/*---------------------------------------------------------------------------*/
	/* If rci_request=anything else, then dcg subroutine failed                  */
	/* to compute the solution vector: solution[n]                               */
	/*---------------------------------------------------------------------------*/
	goto failure;
	/*---------------------------------------------------------------------------*/
	/* Reverse Communication ends here                                           */
	/* Get the current iteration number into itercount                           */
	/*---------------------------------------------------------------------------*/
getsln: dcg_get(&n,solution,rhs,&rci_request,ipar,dpar,tmp,&itercount);
	//fprintf(stderr,"code %d\n", rci_request);
	/*---------------------------------------------------------------------------*/
	/* Print solution vector: solution[n] and number of iterations: itercount    */
	/*---------------------------------------------------------------------------*/
	//printf("The system has been solved\n");
	//printf("\nNumber of iterations: %d\n",itercount);

	/*-------------------------------------------------------------------------*/
	/* Release internal MKL memory that might be used for computations         */
	/* NOTE: It is important to call the routine below to avoid memory leaks   */
	/* unless you disable MKL Memory Manager                                   */
	/*-------------------------------------------------------------------------*/
	//MKL_FreeBuffers();
        mkl_free_buffers();
	free(tmp);

	
    //fprintf(stdout,"This example has successfully PASSED through all steps of computation!\n");
    return itercount;
	/*-------------------------------------------------------------------------*/
	/* Release internal MKL memory that might be used for computations         */
	/* NOTE: It is important to call the routine below to avoid memory leaks   */
	/* unless you disable MKL Memory Manager                                   */
	/*-------------------------------------------------------------------------*/
failure:
        mkl_free_buffers();
	free(tmp);
				 fprintf(stderr,"cg FAILED as the solver has returned the ERROR "); 
				 fprintf(stderr,"code %d", rci_request);
         return - 1;
}

