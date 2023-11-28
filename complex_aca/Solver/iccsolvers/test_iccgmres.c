/********************************************************************************
/   Content:
/   Intel MKL RCI (P)FGMRES ((Preconditioned) Flexible Generalized Minimal
/                                                        RESidual method) example
/ *******************************************************************************/

/*---------------------------------------------------------------------------
/   Example program for solving non-symmetric indefinite system of equations
/   Fully advanced case: full functionality of RCI FGMRES solver is exploited
/ ---------------------------------------------------------------------------*/

#include <stdio.h>
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"



#define N 5
#define size 128

// data structure for matrix-vector product
typedef struct smvData {
	MKL_INT n;
	MKL_INT *ia,*ja;
	double* a;

} mvData;

// matrix-vector multiplication function
void matvec(double *x, double *y, mvData *mvd) {
	double dvar;
	char cvar;
	cvar='N';
	mkl_dcsrgemv(&cvar, &(mvd->n), mvd->a, mvd->ia, mvd->ja, x, y);
}

// printing function for convergence monitoring
void prnfnc(int itn, double err) {
	fprintf(stdout,"iter= %6d err= %le\n",itn,err);
	fflush(stdout);
}

int main(void) {

	/*---------------------------------------------------------------------------
	/  Define arrays for the upper triangle of the coefficient matrix
	/  Compressed sparse row storage is used for sparse representation
	/ ---------------------------------------------------------------------------*/
	MKL_INT ia[6]={1,3,6,9,12,14};
	MKL_INT ja[13]={    1,        3,
						 1,   2,        4,
								2,   3,        5,
									  3,   4,   5,
											 4,   5  };
	double A[13]={ 1.0,     -1.0,
					  -1.0, 1.0,     -1.0,
					        1.0,-2.0,      1.0,
					            -1.0, 2.0,-1.0,
					                 -1.0,-3.0 };
	mvData mvd;
	/*---------------------------------------------------------------------------
	/  Allocate storage for the solution/rhs/residual vectors
	/ ---------------------------------------------------------------------------*/
	double expected_solution[N]={-1.0,1.0,0.0,1.0,-1.0};
	double rhs[N];
	double computed_solution[N];
	/*---------------------------------------------------------------------------
	/  Some additional variables to use with the RCI (P)FGMRES solver
	/ ---------------------------------------------------------------------------*/
	MKL_INT itercount, expected_itercount=4;
	MKL_INT RCI_request, i, ivar;
	double dvar;
	char cvar;

	mvd.n=N;
	mvd.ia=ia;
	mvd.ja=ja;
	mvd.a=A;
	printf("--------------------------------------------------------\n\n");
	/*---------------------------------------------------------------------------
	/  Initialize variables and the right hand side through matrix-vector product
	/ ---------------------------------------------------------------------------*/
	ivar=N;
	cvar='N';
	mkl_dcsrgemv(&cvar, &ivar, A, ia, ja, expected_solution, rhs);
	/*---------------------------------------------------------------------------
	/  Initialize the initial guess
	/ ---------------------------------------------------------------------------*/
	for(i=0;i<N;i++)
	{
		computed_solution[i]=0.0;
	}
	// call the gmres solution procedure
	itercount=solve_system_gmres_mv(N, 
			(void (*)(double *, double *, void *)) matvec, (void*) &mvd,
			NULL,  NULL,
			rhs, computed_solution,
			1.0e-6, 10,prnfnc);
	printf(" The system has been solved \n");
	printf("\n The following solution has been obtained: \n");
	for (i=0;i<N;i++) {
		printf("computed_solution[%d]=",i);
		printf("%e\n",computed_solution[i]);
	}
	printf("\n The expected solution is: \n");
	for (i=0;i<N;i++) {
		printf("expected_solution[%d]=",i);
		printf("%e\n",expected_solution[i]);
	}
	printf("\n Number of iterations: %d\n",itercount);
	mkl_free_buffers();
        return 0;
}

