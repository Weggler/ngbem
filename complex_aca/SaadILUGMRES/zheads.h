#include <complex.h>
/* 
| Header file for Complex Solvers package
| Last Modified: Jun. 11, 2007 DOK
*/


typedef struct zSpaFmt *csptr;
typedef struct zSpaFmt {
/*--------------------------------------------- 
| C-style CSR format - used internally
| for all matrices in CSR format 
|---------------------------------------------*/
   int n;
  /*  length of each row */
   int *nzcount;
  /* pointer-to-pointer to store nonzero entries */
   complex double **ma;
  /* pointer-to-pointer to store column indices  */
   int **ja;
} zSparMat;

/* -------------------------------------------------------------------*/
typedef struct zPerMat4 *p4ptr;
typedef struct zPerMat4 {
/*------------------------------------------------------------
| struct for storing the block LU factorization 
| contains all the block factors except the 
| data related to the last block. 
| n       = size of current block
| symperm = whether or not permutations are symmetric.
|           used only in cleanP4..
| nB      = size of B-block
| L, U    = ILU factors of B-block
| F, E    = sparse matrices in (1,2) and (2,1) 
|           parts of matrix. 
| perm    = (symmetric) permutation used in factorization
|           comes from the independent set ordering
| rperm   = unsymmetric permutation (rows) not used in this
|           version -- but left here for compatibility..
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|----------------------------------------------------------*/ 
  int n;                  
  int nB; 
  int symperm;
  /*   LU factors  */
  struct zSpaFmt *L;
  struct zSpaFmt *U;
  /* E, F blocks   */
  struct zSpaFmt *E;
  struct zSpaFmt *F;
  int *rperm;       /* row permutation    */ 
  int *perm;        /* col. permutation   */
  double *D1 ;
  double *D2 ;
  complex double *wk;
  /* pointer to next and previous struct */
  p4ptr prev; 
  p4ptr next;
} zPer4Mat; 
/* -------------------------------------------------------------------*/
typedef struct zILUTfac *ilutptr;
typedef struct zILUTfac {
/*------------------------------------------------------------
| struct for storing data related to the last schur complement 
| we need to store the C matrix associated with the last block
| and the ILUT factorization of the related Schur complement.
| 
| n       = size of C block = size of Schur complement
| C       = C block of last level matrix. 
| L, U    = ILU factors of last schur complement. 
|
| meth[4] = parameters for defining variants in factorization 
|           - see function readin for details
| rperm    = row permutation used for very nonsymmetric matrices 
|            [such as bottleneck transversal] -- NOT IN THIS VERSION
| perm2     = unsymmetric permutation (columns) - used primarily
|           for the ILUTP version of ILUT/.. 
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|-----------------------------------------------------------*/
   int n;                  
  /*-------------------- C matrix of last block */
   struct zSpaFmt *C;
  /* LU factorization       */
   struct zSpaFmt *L;
   struct zSpaFmt *U;
  /*--------------------  related to variants and methods */
  /*    int meth[4];   */
  int *rperm;   /* row single-sinded permutation */
  int *perm;    /* column perm .                */
  int *perm2;   /* column permutation coming from pivoting in ILU */ 
   double *D1;
   double *D2;
   complex double *wk;
} zIluSpar;

typedef struct zarms_st *arms;
typedef struct zarms_st {
  /* this is the arms preconditioner struct 
  | it consists of a linked list of p4mat structs
  | and the ILUt factorization (in the form of an 
  | IluSpar struct  
  |---------------------------------------------- */
  int n;  
  int nlev;
  ilutptr ilus;
  p4ptr levmat;
} zarmsMat;

typedef int idxtype;

typedef struct Graph_t *GraphPtr; 
typedef struct Graph_t {
  /* simplied version of metis graph == will be updated later to match metis
     definition so as to use nested dissections */
  int nvtxs, nedges;	  /* The # of vertices and edges in the graph */
  idxtype *xadj;	  /* Pointers to the locally stored vertices */
  idxtype *adj;	          /* Array storing the adjacency lists */
  idxtype *label;
} GType; 
/*-------------------- old-style csr matrix struct                      */
typedef struct zCSRregular *csrCmat;
typedef struct zCSRregular{
complex double *a;
int *ja;
int *ia;
} zCSRMat;

typedef struct zILUfac {
    int n;
    csptr L;      /* L part elements                            */
    complex double *D;    /* diagonal elements                  */
    csptr U;      /* U part elements                            */
    int *work;    /* working buffer */
} zILUSpar, zLDUmat, *iluptr;



/*------ Combined matrix and Preconditioner structures -----*/
typedef struct _SMat {
  /*-------------------- 2 types of matrices so far */
  int n; 
  int Mtype;           /*--  type 1 = CSR, 2 = LDU    */
  csptr CSR;           /* place holder for a CSR/CSC type matrix */
  iluptr LDU;          /* struct for an LDU type matrix          */
  void (*zmatvec)(struct _SMat*, complex double *, complex double *);
} SMat, *SMatptr;

typedef struct _SPre {
  /*-------------------- 2 types of matrices so far */
  int Ptype;           /*-- Ptype 1 = ILU, 2 = Crout */
  iluptr   ILU;        /* struct for an ILU type preconditioner */
  arms ARMS;           /* struct for ARMS     preconditioner */
  int (*zprecon) (complex double *, complex double *, struct _SPre*); 
} SPre, *SPreptr;


