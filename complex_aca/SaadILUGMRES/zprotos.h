 /*---------------------------------------------------------------------------* 
 * all prototypes are grouped in this file 
 *---------------------------------------------------------------------------*/
#include <complex.h>
#include "zdefs.h"

void errexit(char *f_str, ...);

/* zsets */
void *Malloc(int nbytes, char *msg );
int zmallocRow( iluptr lu, int nrow );
int zsetupCS(csptr amat, int len);
int zcleanCS(csptr amat);
int znnzCS( csptr amat );
int zcscpy(csptr amat, csptr bmat);
int zsetupILU( iluptr lu, int n );
int zcleanILU( iluptr lu );
int zsetupP4 (p4ptr amat, int Bn, int Cn,  csptr F,  csptr E);
int zcleanP4(p4ptr amat);
int zsetupILUT(ilutptr amat, int len);
int zcleanILUT(ilutptr amat, int indic);
void zsetup_arms (arms Levmat);
int zcleanARMS(arms ArmsPre);
int zcsSplit4(csptr amat, int bsize, int csize, csptr B, csptr F,
	     csptr E, csptr C);
int zCSRcs(int n, complex double *a, int *ja, int *ia, csptr bmat);
int zCOOcs(int n, int nnz, complex double *a, int *ja, int *ia, csptr bmat); 
int zlev4_nnz(p4ptr levmat, int *lev, FILE *ft);
int znnz_arms (arms PrecMat, int nlev, FILE *ft);
int znnz_ilu(iluptr lu, FILE *ft );
int outputLU( iluptr lu, char *filename );
int zCSClum( int n, complex double *a, int *ja, int *ia, iluptr mat, int typ);

/* zMatOps */
void zmatvec(csptr mata, complex double *x, complex double *y);
void zLsol(csptr mata, complex double *b, complex double *x);
void zUsol(csptr mata, complex double *b, complex double *x);
int zdescend(p4ptr levmat, complex double *x, complex double *wk);
int zascend (p4ptr levmat, complex double *x, complex double *wk);
void zmatvecz(csptr mata, complex double *x, complex double *y, complex double *z);
p4ptr zLvsol2(complex double *x, int nlev, p4ptr levmat, ilutptr ilusch);
int zUvsol2(complex double *x, int nlev, int n, p4ptr levmat,
	   ilutptr ilusch);
int zarmsol2(complex double *x,  arms Prec);
void zSchLsol(ilutptr ilusch, complex double *y);
void zSchUsol(ilutptr ilusch, complex double *y);
int zlusolC( complex double *y, complex double *x, iluptr lu );
int zlumsolC(complex double *y, complex double *x, iluptr lu );
int zrpermC(csptr mat, int *perm);
int zcpermC(csptr mat, int *perm);
int zdpermC(csptr mat, int *perm);
void zlumatvec( iluptr mat,complex double *x,complex double *y );
int zcondestLU( iluptr lu, complex double *y, complex double *x, FILE *fp );
int zcondestArms(arms armspre, complex double *y, FILE *fp );
void zmatvecCSR(SMatptr mat, complex double *x, complex double *y);
void zmatvecLDU(SMatptr mat, complex double *x, complex double *y);
int  zpreconILU(complex double *x, complex double *y, SPreptr mat);
int  zpreconLDU(complex double *x, complex double *y, SPreptr mat);
int  zpreconARMS(complex double *x, complex double *y, SPreptr mat);
int zfgmres(SMatptr Amat, SPreptr PreMat, complex double *rhs, 
	    complex double *sol, double tol, int im, int *itmax, FILE *fp);

/* zmisc */ 
int zqsplitC(complex double *a, int *ind, int n, int ncut);
void qsplitC(double *a, int *ind, int n, int Ncut);
int zSparTran(csptr amat, csptr bmat, int job, int flag);
void zswapj(int v[], int i, int j);
void zswapm(complex double v[], int i, int j);
void swapm(double v[], int i, int j);
void zqsortC(int *ja, complex double *ma, int left, int right, int abval);
int zroscalC(csptr mata, double *diag, int nrm);
int zcoscalC(csptr mata, double *diag, int nrm);
void zdscale(int n, double *dd, complex double *x, complex double * y);
void zprintmat(FILE *ft, csptr A, int i0, int i1);
void qsortR2I(double *wa, int *cor1, int *cor2, int left, int right);
void zqsort2C(int *ja, complex double *ma, int left, int right, int abval);
void zqqsort(int *ja, complex double *ma, int left, int right);
void zhilosort(csptr mat, int abval, int hilo);
void zqsort3i(int *wa, int *cor1, int *cor2, int left, int right);

/* zindsetC */
int zadd2is(int *last, int nod, int *iord, int *riord); 
int zadd2com(int *nback, int nod, int *iord, int *riord); 
int zindsetC(csptr mat, int bsize, int *iord, int *nnod, double tol); 
int zweightsC(csptr mat, double *w);

/* zilutpC */         
int zilutpC(csptr amat, double *droptol, int *lfil, double permtol,
		int mband, ilutptr ilusch);
int zilutD(csptr amat, double *droptol, int *lfil, ilutptr ilusch);

/* zpiluNEW */
int zpilu(p4ptr amat, csptr B, csptr C, double *droptol, 
	 int *lfil, csptr schur);

/* zPQ */
int zPQperm(csptr mat, int bsize, int *Pord, int *Qord, int *nnod, 
	     double tol);
int zpreSel(csptr mat, int *icor, int *jcor, int job, double tol, int *count);

/* zprtC */
int zprtC(csptr Amat, int io);

/* precons */
int zarms2(csptr Amat, int *ipar, double *droptol, int *lfil, 
	  double tolind, arms PreMat, FILE *ft);

int zilukC( int lofM, csptr csmat, iluptr lu, FILE *fp );

int zilut(csptr csmat, iluptr lu, int lfil, double tol, FILE *fp );

int zilutc(iluptr mt, iluptr lu, int lfil, double tol, int drop_type, FILE *fp );

/* auxill */
void zoutput_perm( int n, int *perm, FILE *f );

/* ztools */
void zroscal (int* nrow, int* job, int* nrm, complex double *a, int *ja, 
	     int *ia, double *diag, complex double *b, int *jb, int *ib, int *ierr) ;
void zcoscal (int* nrow, int* job, int* nrm, complex double *a, int *ja, 
	     int *ia, double *diag, complex double *b, int *jb, int *ib, int *ierr) ;
void zreadmtc(int*,  int*,  int*,  char*,  complex double*,  int*,  int*,  complex double*, 
	     int*,  char*,  int*,  int*,  int*,  char*,  char*, 
	     char*,  int*) ;
void zcsrcsc (int*, int*, int*, complex double*, int*, int*, complex double*, int*, int*) ;
void qsplit(double *a, int *ind, int *n, int *Ncut);
double dznrm2(int *, complex double *, int *);
//complex double zdotc(int *, complex double *, int *, complex double *, int *);
void zaxpy(int *, complex double *, complex double *, int *, complex double *, int *);

/* givens */
void zclartg(complex double f, complex double g, double *cs, complex double *sn, complex double *rot);



/* end protos */ 
