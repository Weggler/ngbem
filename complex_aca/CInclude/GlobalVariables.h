//
// S. Rjasanow: (Matrix) Adaptive Cross Approximation
//
//------------------------------------------------------------------
//
//  header name   -  GlobalVariables
//
//------------------------------------------------------------------
//
//  last revision -  Mar 12
//  purpose       -  global variables                               
//
//------------------------------------------------------------------

# include "Includes.h"

// NCluMin    - min number of cluster 
// ACAEta     -
// ACABeta    -
// ACAGamma   -
// InfoLevel  -
long NCluMin;
double ACAEta;
double ACAAlpha;
double ACABeta;
double ACAGamma;
long InfoLevel;

// SymmType     - symmetry-type flag
// MaxRank      - max rank of low rank approximation
// ACAEps       - minimal accuracy of approxmation
// Memory       - 
// AppMemory    - 
SymmType Symmetry;
long MaxRank;
double ACAEps;
double Memory,AppMemory;

// PrecondType  - precondition-type flag
// NMaxILU      - something with incomplete LU
// LFillILU     - something with fill in when incomplete LU
// LevelBE      - connectivity level of basis elements (Precond==BasisElement)
// LevelNB      - connectivity level of elements (Precond==BasisElement)
PrecondType Precond;
long NMaxILU,LFillILU;
long LevelBE;
long LevelNB;

// HMax,HMin     - max,min mesh width
// CB            - ?
// Count1,Count2 - counting indices, help variables
double HMax,HMin,CB;
long Count1, Count2;

/* Environment */
char Prefix[4096];
char Suffix[4096];
char User[4096];

/* Global Hierarchical Block Matrix */
// BockHMatrix   - global hierarchical block matrix
// ACAParameters - set ACA parameters (specifying the approximation quality...)
// NCLUSTERS     - cluster flags (the flag is a number ?) 
// CLUSTERS      - list of sets of NCLUSTER clusters (for instance for vertices, medge, mdle) 
//                 here, I use two (Maxwell): a cluster wrt medge nodes
//                 and one wrt mdle nodes 
// PERMU         - list of permutation sets for the clusters
// PAIRS         - list of sets of corresponding  cluster pairs
// NMAX,NMAX     - maximal row, column?
//double complex *EBlock_new;
BlockHMatrix *BHMAT;
ACAParameters PARAMACA;
long NCLUSTERS[2];
Cluster *CLUSTERS[2];
long *PERMU[2];
Pair *PAIRS[4];
long NPAIRS[4];
long NMAX[4],MMAX[4];

/* more for MACA */
// NLoc          - the element matrix dimensions for 2 Clusters
long *NLOC[2];
int *NR_CON;
int *NODES_CON;
int NDOFTOT;
int NDOFMACA;
