#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <complex.h>
#include <assert.h>
#include "zheads.h"

void errexit( char *f_str, ... ){
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s\n", out1);

  fprintf(stdout, out2);
  fflush(stdout);

  exit( -1 );
}

void *Malloc(int nbytes, char *msg )
{ /* allocates space -- or exits in case of
     failure */
  void *ptr;
  if (nbytes == 0)
    return NULL;
  ptr = (void *) malloc(nbytes);
  if (ptr == NULL) {
    fprintf(stderr,"Mem. alloc. ERROR in %s. Requested bytes: %d bytes",
	    msg, nbytes );
     exit( -1 );
  }
  return ptr;
}

int zmallocRow( iluptr lu, int nrow )
{
/*----------------------------------------------------------------------
| Prepare space of a row according to the result of level structure
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|     nrow  =  the current row to deal with
|
| On return:
|===========
|
|    lu->L->ma[nrow][...]
|      ->U->ma[nrow][...]
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int nzcount = lu->L->nzcount[nrow];
    lu->L->ma[nrow] = (complex double *)Malloc( sizeof(complex double)*nzcount, "mallocRow" );
    nzcount = lu->U->nzcount[nrow];
    lu->U->ma[nrow] = (complex double *)Malloc( sizeof(complex double)*nzcount, "mallocRow" );
    return 0;
}
/*---------------------------------------------------------------------
|     end of mallocRow
|--------------------------------------------------------------------*/


int zsetupCS(csptr amat, int len)
{
/*---------------------------------------------------------------------- 
| Initialize SpaFmt structs.
|----------------------------------------------------------------------
| on entry: 
|========== 
| ( amat )  =  Pointer to a SpaFmt struct.
|     len   =  size of matrix
|
| On return:
|===========
|
|  amat->n
|      ->*nzcount
|      ->**ja
|      ->**ma
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   amat->n = len;
   amat->nzcount = (int *) Malloc(len*sizeof(int), "setupCS:1" );
   amat->ja = (int **) Malloc(len*sizeof(int *), "setupCS:2" );
   amat->ma = (complex double **) Malloc(len*sizeof(complex double *), "setupCS:3" );
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupCS
|--------------------------------------------------------------------*/

int zcleanCS(csptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for SpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SpaFmt struct.
|     len   =  size of matrix
|--------------------------------------------------------------------*/
   /*	*/
  int i; 
  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;

  for (i=0; i<amat->n; i++) {
    if (amat->nzcount[i] > 0) {
      if (amat->ma[i]) free(amat->ma[i]);
      if (amat->ja[i]) free(amat->ja[i]);
    }
  }	
  if (amat->ma) {
    free(amat->ma);
    amat->ma = NULL;
  }
  if (amat->ja) {
    free(amat->ja);
    amat->ja = NULL;
  }
  if (amat->nzcount) {
    free(amat->nzcount); 
    amat->nzcount = NULL; 
  }
  if (amat){
    free(amat);  
    amat = NULL;
  }
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanCS
|--------------------------------------------------------------------*/

int znnzCS( csptr amat )
{
/* counts the number of nonzero elements in CSR matrix A */
    int nnz = 0, i, n = amat->n;
    for( i = 0; i < n; i++ ) {
        nnz += amat->nzcount[i];
    }
    return nnz;
}

int zcscpy(csptr amat, csptr bmat)
{
/*----------------------------------------------------------------------
| Convert CSR matrix to SpaFmt struct
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )   = Matrix stored in SpaFmt format
|
|
| On return:
|===========
|
| ( bmat )  =  Matrix stored as SpaFmt struct containing a copy
|              of amat 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int j, len, size=amat->n;
   complex double *bma;
   int *bja;
/*------------------------------------------------------------*/
   for (j=0; j<size; j++) {
     len = bmat->nzcount[j] = amat->nzcount[j];
     if (len > 0) {
       bja = (int *) Malloc(len*sizeof(int), "cscpy:1" );
       bma = (complex double *) Malloc(len*sizeof(complex double), "cscpy:2" );
       memcpy(bja,amat->ja[j],len*sizeof(int));
       memcpy(bma,amat->ma[j],len*sizeof(complex double));
       bmat->ja[j] = bja;
       bmat->ma[j] = bma;
      }	
   }
   return 0;
}
/*-----------------------------------------------------------------------*/

int zsetupILU( iluptr lu, int n )
{
/*----------------------------------------------------------------------
| Initialize ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|       n   =  size of matrix
|
| On return:
|===========
|
|    lu->n
|      ->L     L matrix, SpaFmt format
|      ->D     Diagonals
|      ->U     U matrix, SpaFmt format
|      ->work  working buffer of length n
|      ->bf    buffer
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    lu->n  = n;
    lu->D = (complex double *)Malloc( sizeof(complex double) * n, "setupILU" );
    lu->L = (csptr)Malloc( sizeof(zSparMat), "setupILU" );
    zsetupCS( lu->L, n);
    lu->U = (csptr)Malloc( sizeof(zSparMat), "setupILU" );
    zsetupCS( lu->U, n);
    lu->work = (int *)Malloc( sizeof(int) * n, "setupILU" );
    return 0;
}
/*---------------------------------------------------------------------
|     end of setupILU
|--------------------------------------------------------------------*/


int zcleanILU( iluptr lu )
{
/*----------------------------------------------------------------------
| Free up memory allocated for ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|--------------------------------------------------------------------*/
  if( NULL == lu ) return 0;
  if( lu->D ) {
    free( lu->D );
  }
  zcleanCS( lu->L );
  zcleanCS( lu->U );  
  if( lu->work ) free( lu->work );
  free( lu );
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanILU
|--------------------------------------------------------------------*/ 
int zsetupP4 (p4ptr amat, int Bn, int Cn,  csptr F,  csptr E) 
{
/*----------------------------------------------------------------------
| initialize PerMat4 struct given the F, E, blocks.  
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a PerMat4 struct.
|     Bn    =  size of B block
|     Cn    =  size of C block
|     F, E  = the two blocks to be assigned to srtruct - without the
|
| On return:
|===========
|
|  amat->L                for each block: amat->M->n
|      ->U                                       ->nzcount
|      ->E                                       ->ja
|      ->F                                       ->ma
|      ->perm
|      ->rperm       (if meth[1] > 0)
|      ->D1          (if meth[2] > 0)
|      ->D2          (if meth[3] > 0)
|
|  Scaling arrays are initialized to 1.0.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int n;
   int zsetupCS(csptr, int);
   /* size n */
   n = amat->n = Bn + Cn;
   amat->nB = Bn; 
/* amat->perm = (int *) Malloc(n*sizeof(int), "setupP4:1" ); */
/*   assign space for wk -- note that this is only done at 1st level
     at other levels, copy pointer of wk from previous level */
   if (amat->prev == NULL)  /* wk has 2 * n entries now */
     amat->wk   = (complex double *) Malloc(2*n*sizeof(complex double), "setupP4:2" );
   else 
     amat->wk = (amat->prev)->wk; 

/*-------------------- L and U */ 
   amat->L = (csptr) Malloc(sizeof(zSparMat), "setupP4:3" );
   if (zsetupCS(amat->L, Bn)) return 1;
   /*    fprintf(stdout,"  -- BN %d   Cn   %d \n", Bn,Cn);  */
   amat->U = (csptr) Malloc(sizeof(zSparMat), "setupP4:4" );
   if (zsetupCS(amat->U, Bn)) return 1;

   amat->F = F; 
   amat->E = E; 
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupP4 
|--------------------------------------------------------------------*/

int zcleanP4(p4ptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for Per4Mat structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a Per4Mat struct.
|--------------------------------------------------------------------*/
  int zcleanCS(csptr);
/*  -------------------------- */
  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;
   

  if (amat->perm) {
    if (amat->perm) free(amat->perm); 
    amat->perm = NULL;
  }
  
  if (!amat->symperm) { 
    if (amat->rperm) free(amat->rperm); 
    amat->rperm = NULL;
  } 
  
  if (amat->F) {
    zcleanCS(amat->F); 
    amat->F = NULL;
  }
  if (amat->E) {
    zcleanCS(amat->E); 
    amat->E = NULL;
  }
  if (amat->L) {
    zcleanCS(amat->L);
    amat->L = NULL;
   }
  if (amat->U) {
    zcleanCS(amat->U);
    amat->U = NULL;
  }
  
  if (amat->prev == NULL) 
    if (amat->wk) free(amat->wk);  
  
  if (amat->D1) free(amat->D1);
  if (amat->D2) free(amat->D2);
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanP4
|--------------------------------------------------------------------*/

int zsetupILUT(ilutptr amat, int len)
{
/*----------------------------------------------------------------------
| Allocate pointers for ILUTfac structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a ILUTfac struct.
|     len   =  size of L U  blocks
|
| On return:
|===========
|
|  amat->L                for each block: amat->M->n
|      ->U                                       ->nzcount
|                                                ->ja
|                                                ->ma
|      ->rperm       (if meth[0] > 0)
|      ->perm2       (if meth[1] > 0)
|      ->D1          (if meth[2] > 0)
|      ->D2          (if meth[3] > 0)
|
|  Permutation arrays are initialized to the identity.
|  Scaling arrays are initialized to 1.0.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int zsetupCS(csptr, int);
  amat->n = len;
 amat->wk = (complex double *) Malloc(2*len*sizeof(complex double), "setupILUT:5" );
  amat->L = (csptr) Malloc(sizeof(zSparMat), "setupILUT:6" );
  if (zsetupCS(amat->L, len)) return 1;
  amat->U = (csptr) Malloc(sizeof(zSparMat), "setupILUT:7" );
  if (zsetupCS(amat->U, len)) return 1;
  return 0;
}    
/*---------------------------------------------------------------------
|     end of setupILUT
|--------------------------------------------------------------------*/
int zcleanILUT(ilutptr amat, int indic)
{
/*----------------------------------------------------------------------
| Free up memory allocated for IluSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a IluSpar struct.
|  indic    = indicator for number of levels.  indic=0 -> zero level.
|--------------------------------------------------------------------*/
  int zcleanCS(csptr);
  /*----------------*/
   
  if (amat->wk) {
    free(amat->wk); 
    amat->wk = NULL;
  }
  zcleanCS(amat->L);
  zcleanCS(amat->U);

  if (indic) zcleanCS(amat->C);  
/*-------------------- nonsymmetric permutation */
  if (amat->rperm) {
    free(amat->rperm);
    amat->rperm = NULL;
  }
  if (amat->perm) {
    free(amat->perm); 
    amat->perm = NULL;
  }
  
/*-------------------- ilutp permutation */
  if (amat->perm2) free(amat->perm2);
/*-------------------- diagonal scalings */
   if (amat->D1) free(amat->D1);
   if (amat->D2) free(amat->D2);
   return 0;
}
/*---------------------------------------------------------------------
|     end of cleanILUT
|--------------------------------------------------------------------*/


void zsetup_arms (arms Levmat) {
  Levmat->ilus = (ilutptr) Malloc(sizeof(zIluSpar), "setup_arms:ilus" );
  Levmat->levmat = (p4ptr) Malloc(sizeof(zPer4Mat), "setup_arms:levmat" );
}
int zcleanARMS(arms ArmsPre)
{
  p4ptr amat = ArmsPre->levmat;
  ilutptr cmat = ArmsPre->ilus;
/*----------------------------------------------------------------------
| Free up memory allocated for entire ARMS preconditioner.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a Per4Mat struct.
| ( cmat )  =  Pointer to a IluSpar struct.
|--------------------------------------------------------------------*/
/* case when nlev == 0 */  
  int indic=(amat->nB != 0) ;
    /*  && amat->next !=NULL) ; */
  
  p4ptr levc, levn;
  int zcleanCS(csptr);
  int zcleanP4(p4ptr);
  int zcleanILUT(ilutptr, int);

  levc = amat; 

  if (indic) { 
    while (levc) {
      if (zcleanP4(levc)) return(1) ; 
      levn = levc->next;
      free(levc);
      levc = levn;
    }		
  }	
   else 	
     if (amat) {
       free(amat) ; 
       amat = NULL;
     }
  
  zcleanILUT(cmat,indic); 
  
  
  if (cmat) {
    free(cmat);	
    cmat = NULL;
  }
  
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanARMS 
|--------------------------------------------------------------------*/

int zcsSplit4(csptr amat, int bsize, int csize, csptr B, csptr F,
	     csptr E, csptr C)
{
/*---------------------------------------------------------------------
| Convert permuted csrmat struct to PerMat4 struct 
|                - matrix already permuted
|----------------------------------------------------------------------
| on entry:
|========== 
| ( amat )  =  Matrix stored in SpaFmt format.
|              Internal pointers (and associated memory) destroyed before
|              return.
|
| On return:
|===========
|
| B, E, F, C = 4 blocks in 
| 
|          | B   F |      
|   Amat = |       | 
|          | E   C | 
| 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int j, j1, numr, numl, ind, newj, rowz, *rowj, *new1j, *new2j;
   complex double *rowm, *new1m, *new2m;
/*---------------------------------------------------------------------
|     Sort the matrix and separate into   |  B  F  |
|                                         |        |
|                                         |  E  C  |
|--------------------------------------------------------------------*/
   if (zsetupCS(B,bsize)) goto label111; 
   if (zsetupCS(F,bsize)) goto label111;
   if (zsetupCS(E,csize)) goto label111;
   if (zsetupCS(C,csize)) goto label111;
   new1j = (int *) Malloc(bsize*sizeof(int), "csSplit4:1" );
   new2j = (int *) Malloc(csize*sizeof(int), "csSplit4:2" );
   new1m = (complex double *) Malloc(bsize*sizeof(complex double), "csSplit4:3" );
   new2m = (complex double *) Malloc(csize*sizeof(complex double), "csSplit4:4" );
/*    B and F blocks */ 
   for (j=0; j<bsize; j++) {
      numl = numr = 0;
      rowz = amat->nzcount[j];
      rowj = amat->ja[j];
      rowm = amat->ma[j];
      for (j1=0; j1<rowz; j1++) {
	 if (rowj[j1]<bsize) numl++;
	 else numr++;
      }
      B->nzcount[j] = numl;
      F->nzcount[j] = numr;
      if (numl>0) {
	 B->ja[j] = (int *) Malloc(numl*sizeof(int), "csSplit4:5" );
	 B->ma[j] = (complex double *) Malloc(numl*sizeof(complex double), "csSplit4:6" );
      }
      if (numr>0) {
	 F->ja[j] = (int *) Malloc(numr*sizeof(int), "csSplit4:7" );
	 F->ma[j] = (complex double *) Malloc(numr*sizeof(complex double), "csSplit4:8" );
      }
      numl = numr = 0;
      for (j1=0; j1<rowz; j1++) {
	 newj = rowj[j1];
	 if (newj<bsize) {
	    new1j[numl] = newj;
	    new1m[numl] = rowm[j1];
	    numl++;
	 }
	 else {
	    new2j[numr] = newj - bsize;
	    new2m[numr] = rowm[j1];
	    numr++;
	 }
      }
      memcpy(B->ja[j], new1j, numl*sizeof(int));
      memcpy(B->ma[j], new1m, numl*sizeof(complex double));
      memcpy(F->ja[j], new2j, numr*sizeof(int));
      memcpy(F->ma[j], new2m, numr*sizeof(complex double));
   }
/*    E and C blocks */
   for (j=0; j<csize; j++) {
      numl = numr = 0;
      ind = bsize + j;
      rowz = amat->nzcount[ind];
      rowj = amat->ja[ind];
      rowm = amat->ma[ind];
      for (j1=0; j1<rowz; j1++) {
	 if (rowj[j1]<bsize) numl++;
	 else numr++;
      }
      E->nzcount[j] = numl;
      C->nzcount[j] = numr;
      if (numl>0) {
	E->ja[j] = (int *) Malloc(numl*sizeof(int), "csSplit4:9" );
	E->ma[j] = (complex double *) Malloc(numl*sizeof(complex double), "csSplit4:10" );
      }	
      if (numr>0) {
	C->ja[j] = (int *) Malloc(numr*sizeof(int), "csSplit4:11" );
	C->ma[j] = (complex double *) Malloc(numr*sizeof(complex double), "csSplit4:12" );
      }		
      numl = numr = 0;
      for (j1=0; j1<rowz; j1++) {
	newj = rowj[j1];
	if (newj<bsize) {
	  new1j[numl] = newj;
	  new1m[numl] = rowm[j1];
	  numl++;
	}
	else {
	  new2j[numr] = newj - bsize;
	  new2m[numr] = rowm[j1];
	  numr++;
	}
      }
      memcpy(E->ja[j], new1j, numl*sizeof(int));
      memcpy(E->ma[j], new1m, numl*sizeof(complex double));
      memcpy(C->ja[j], new2j, numr*sizeof(int));
      memcpy(C->ma[j], new2m, numr*sizeof(complex double));
   }

   if (new1j) free(new1j);
   if (new2j) free(new2j);
   if (new1m) free(new1m);
   if (new2m) free(new2m);
   return 0;
label111:
   return 1;
}
/*---------------------------------------------------------------------
|     end of csSplit4
|--------------------------------------------------------------------*/

int zCSRcs(int n, complex double *a, int *ja, int *ia, csptr bmat)
{
/*----------------------------------------------------------------------
| Convert CSR matrix to SpaFmt struct
|----------------------------------------------------------------------
| on entry:
|==========
| a, ja, ia  = Matrix stored in CSR format (with FORTRAN indexing).
|
| On return:
|===========
|
| ( bmat )  =  Matrix stored as SpaFmt struct.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int i, j, j1, len, st; 
  complex double *bra;
  int *bja;
  int zsetupCS(csptr, int);
  /*	setup data structure for bmat (csptr) struct */
  if (zsetupCS(bmat, n)) {
    printf(" ERROR SETTING UP bmat IN SETUPCS \n") ;
    exit(0);
  }
  /*				*/
  st = ia[0];
  for (j=0; j<n; j++) {
    len = ia[j+1] - ia[j];
    bmat->nzcount[j] = len;
    if (len > 0) {
      bja = (int *) Malloc(len*sizeof(int), "CSRcs:1" );
      bra = (complex double *) Malloc(len*sizeof(complex double), "CSRcs:2" );
      i = 0;
      for (j1=ia[j]-st; j1<ia[j+1]-st; j1++) {
	bja[i] = ja[j1] - st;
	bra[i] = a[j1] ;  
	i++;
      }
      bmat->ja[j] = bja;
      bmat->ma[j] = bra;
    }
  }	
  return 0;
}
/*---------------------------------------------------------------------
|     end of CSRcs
|--------------------------------------------------------------------*/

/* %%% */

int zCOOcs(int n, int nnz, complex double *a, int *ja, int *ia, csptr bmat)
{
/*----------------------------------------------------------------------
| Convert COO matrix to SpaFmt struct
|----------------------------------------------------------------------
| on entry:
|==========
| a, ja, ia  = Matrix stored in COO format  -- a = complex entries
|                                             ja = column indices
|                                             ia = row indices 
| On return:
|===========
|
| ( bmat )  =  Matrix stored as SpaFmt struct.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int i,  k, k1, l;
  int *len;
  int zsetupCS(csptr, int);
  /*-------------------- setup data structure for bmat (csptr) struct */
  if (zsetupCS(bmat, n)) {
    printf(" ERROR SETTING UP bmat IN SETUPCS \n") ;
    exit(0);
  }
 /*				*/
/*-------------------- determine lengths */
   len =  (int*) Malloc(n*sizeof(int), "COOcs:0" );
   assert(len != NULL);

   for (k=0; k<n; k++) 
     len[k] = 0;
   for (k=0; k<nnz; k++)
      ++len[ia[k]]; 
/*-------------------- allocate          */
   for (k=0; k<n; k++) {
     l = len[k];
     bmat->nzcount[k] = l;
     if (l > 0) {
       bmat->ja[k] = (int *) Malloc(l*sizeof(int), "COOcs:1" );
       bmat->ma[k] = (complex double *) Malloc(l*sizeof(complex double),"COOcs:2" );
     }
     len[k] = 0;
   }
/*-------------------- Fill actual entries */
   for (k=0; k<nnz; k++) {
     i  = ia[k];
     k1 = len[i];
     (bmat->ja[i])[k1] = ja[k] ; 
     (bmat->ma[i])[k1] =  a[k] ;  
     len[i]++; 
   }
   free(len);
   return 0;
}
/*---------------------------------------------------------------------
|     end of COOcs
|--------------------------------------------------------------------*/

/* %%% */



int zlev4_nnz(p4ptr levmat, int *lev, FILE *ft) 
{
  /* counts all nonzero elements in levmat struct  -- 
     recursive */
  int n, nB, nnzT, nnzL, nnzU, nnzF, nnzE, nnzDown=0;
  p4ptr nextmat; 
  n    = levmat->n;
  nB   = (levmat->L)->n; 
  nnzL = znnzCS(levmat->L); 
  nnzU = znnzCS(levmat->U); 
  nnzF = znnzCS(levmat->F); 
  nnzE = znnzCS(levmat->E); 
  nnzT = nnzL+nnzU+nnzF+nnzE;
  /* print */ 
  if (*lev == 0) 
    fprintf(ft,
"\nLev      n     nB    nnzL    nnzU    nnzF    nnzE   subtot\n");  
  fprintf(ft,"%3d %6d %6d %7d %7d %7d %7d %8d\n",
   	  *lev, n, nB, nnzL, nnzU, nnzF, nnzE, nnzT);
  (*lev)++;
  nextmat = levmat->next; 
  if (nextmat != NULL) 
   nnzDown = zlev4_nnz(nextmat, lev, ft);
   return (nnzT+nnzDown); 
}

int znnz_arms (arms PrecMat, int nlev, FILE *ft)
{ 
/*-------------------------------------------------------
| computes and prints out total number of nonzero elements
| used in ARMS factorization 
+--------------------------------------------------------*/
  p4ptr levmat=PrecMat->levmat; ilutptr ilschu = PrecMat->ilus; 
  int ilev=0,nnz_lev,nnz_sch,nnz_tot; 
  nnz_lev = 0; 
  if (nlev) nnz_lev+= zlev4_nnz(levmat, &ilev, ft);
  nnz_sch = znnzCS(ilschu->L)+znnzCS(ilschu->U);
  if (nlev) nnz_sch += znnzCS(ilschu->C);
  nnz_tot = nnz_lev+nnz_sch; 
  fprintf(ft,"\n"); 
  fprintf(ft,"Total nonzeros for interm. blocks.... =  %10d\n",nnz_lev);
  fprintf(ft,"Total nonzeros for last level ....... =  %10d\n",nnz_sch);
  fprintf(ft,"Grand total.......................... =  %10d\n",nnz_tot);
  fprintf(ft,"Size of last Schur complement matrix. =  %10d\n",
	  ilschu->n);
  return nnz_tot;
}

/*-----------------------------------------------------------------------*/
int znnz_ilu(iluptr lu, FILE *ft ){
/*-----------------------------------------------------------------------* 
 * counts the number of nonzero elements in L U factorization            *  
 *-----------------------------------------------------------------------*/
  int n = lu->n, nnz_tot = 0, i;
  int nnz_L=0, nnz_U=0;
  int *nnzcol = lu->L->nzcount, *nzcount = lu->U->nzcount;
  for(i = 0; i < n; i++ ) {
    nnz_L += nnzcol[i];
    nnz_U += nzcount[i];
  }
  nnz_tot = nnz_L + nnz_U + n; 
  if (ft != NULL){
    fprintf(ft,"\n"); 
    fprintf(ft,"Total nonzeros for L block.... =  %10d\n",nnz_L);
    fprintf(ft,"Total nonzeros for D block ... =  %10d\n",n);
    fprintf(ft,"Total nonzeros for U block ... =  %10d\n",nnz_U);
    fprintf(ft,"Grand total................... =  %10d\n",nnz_tot);
  }
  return nnz_tot;
}

int outputLU( iluptr lu, char *filename )
{
/*----------------------------------------------------------------------
| Output the pattern of L\U, [mainly for plotting pattern with matlab] 
----------------------------------------------------------------------*/
    FILE *fmatlab = fopen( filename, "w" );
    int n = lu->n, i, j, nzcount;
    csptr L = lu->L;
    csptr U = lu->U;

    if( !fmatlab ) return -1;
  fprintf( fmatlab, "%d %d 0\n", n, n ); 
  for( i = 0; i < n; i++ ) {
    nzcount = L->nzcount[i];
    for( j = 0; j < nzcount; j++ )
      fprintf( fmatlab, "%d %d 1\n", i+1, L->ja[i][j]+1 );
  }
  for( i = 0; i < n; i++ ) {
    nzcount = U->nzcount[i]; 
    for( j = 0; j < nzcount; j++ )
      fprintf( fmatlab, "%d %d 1\n", i+1, U->ja[i][j]+1 );
  }
  for( i = 0; i < n; i++ )
    fprintf( fmatlab, "%d %d 1\n", i+1, i+1 );
  fclose( fmatlab );
  return 0; 
}

  
