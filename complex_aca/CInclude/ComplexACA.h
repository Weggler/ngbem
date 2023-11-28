//------------------------------------------------------------------
//
//  header name  -  ComplexACA with additional attributes for MACA
//
//------------------------------------------------------------------
//
//  last revision -  Aug 12
//
//   purpose      - definition of C-structures, for instance,
//                    struct Cluster
//                    struct Pair
//                    struct sACAParameters
//                    struct CBlock
//                    struct sHMatrix
//                    struct sBlockHMatrix
//
//------------------------------------------------------------------

/* Saad types -- gmres solver*/
#include "../SaadILUGMRES/zheads.h"

// Level          - cluster level
// Father,Son1,2  - to reconstruct the tree
// Number         - number of points, 
//                  i.e., the dimension of the Cluster 
// PermuPos       - starting position in global PERMU vector to find the 
//                  permuation in this cluster (for a list of clusters)
// Radius         - geometrical
// EVal[3]        - data
// EVec[9]        - for 
// XMin[3]        - its 
// XMax[3]        - construction (?)
// Centre[3]      - from covariance ...
// DiagLength     -
typedef struct sCluster
{
    long Level;
    long Father,Son1,Son2;
    long Number;
    long PermuPos;
    long BasisElement;

    double Radius;
    double EVal[3],EVec[9];
    double XMin[3],XMax[3];
    double Centre[3];
    double DiagLength;
} Cluster;

//  Clu1,2    - the cluster numbers of the pair
//  Type      - its type, i.e.,
//  A11,..,22 - numbers of the CBlocks of the pair 
typedef enum {Admissible=0,Dense=1,Hierarchical=2,Null=3} BlockType;
typedef struct sPair
{
  long Clu1,Clu2;
  BlockType Type;
  long A11,A21,A12,A22;
} Pair;

// NCluMin  - minimal allowed number of cluster points
// MaxRank  - maximal rank of low rank approximation
// NMax,NMax- maximal dimension of ?
// ACAEta   - 
// ACAAlpha -
// ACABeta  -
// ACAGamma -
// ACAKappa -
// ACAEps   - accuracy of approximation
// AppMemory-
// PCMemory -
// Memory   -
// Symmetry - symmetry flag
// Infor    - info flag
typedef enum {NonSymmetric=0,ComplexSymmetric=1} SymmType;
typedef struct sACAParameters
{
  long NCluMin,MaxRank,NMax,MMax;
  double ACAEta,ACAAlpha,ACABeta,ACAGamma,ACAKappa,ACAEps,AppMemory,PCMemory,Memory;
  SymmType Symmetry;
  long Info;
} ACAParameters;

// IBlock,JBlock   - position cluster point permutation vector
// NBlock,MBlock   - dimension of block
// Type            - type of block, see above
// A11,A21,A12,A22 - list of CBlocks non-empty in case
//                   Type == 2
// Data            - ?
// RankClu         - rank wrt to cluster 
// Rank            - rank of the matrices U V^*
// UData,VData     - rectangular matrixes, i.e., low rank approximation 
typedef struct sCBlock
{
  long IBlock,JBlock;
  long NBlock,MBlock;
  BlockType Type;

  struct sCBlock *A11,*A21,*A12,*A22;

  double complex *Data;
  long RankClu;
  long Rank;
  double complex *UData,*VData;
} CBlock;

typedef struct sHBlockInfo
{
  long IBlock,JBlock;
  long NBlock,MBlock;
  BlockType Type;

  long NBlocks,HBlocks,DBlocks,ABlocks;

  double NMemory,HMemory,DMemory,AMemory,Memory;
} HBlockInfo;

typedef struct sSparseCOOMatrix
{
  long NRow,NColumn;
  long NumNZ;

  double complex *Values;
  long *Rows,*Columns;
} SparseCOOMatrix;

// NRow, NColum            - dimension of HMatrix
// NLocRow, NLocColumn     - number of dofs attributed to cluster points
//                           (this is additionally NEEDED FOR MACA)
// NPermuRow, NPermuColumn - total number of cluster points
//                           (this is additionally NEEDED FOR MACA)
// PermuRow, PermuColumn   - permutation of cluster points
// Number                  - ?
// NBlocks                 - Number of CBlocks
// CBlocks                 - list of CBlocks
// Symmetry                - Symmetry flag
// Precond                 - Precondition flag
//                           rest is data specifying the precondition
typedef enum {None=0,Jacobi=1,TriDiag=2,ILU=3,BasisElement=4} PrecondType;
typedef struct sHMatrix
{
  long NRow,NColumn;
  long *NLocRow,*NLocColumn;
  long NPermuRow,NPermuColumn;
  long *PermuRow,*PermuColumn;
  long Number;

  long NBlocks;
  CBlock *CBlocks;

  SymmType Symmetry;
  PrecondType Precond;
  long *PrecondIData;
  double complex *PrecondDData; 
  csptr csmat;
  iluptr lu;
} HMatrix;

typedef struct sGCBlock
{
  long NBlock,MBlock;
  long Cls;
  long *Rows,*Cols;

  BlockType Type;

  double complex *Data;

  long Rank;
  double complex *UData,*VData;
} GCBlock;

typedef struct sGHMatrix
{
  long NRow,NColumn;
  long NMax,MMax;

  long NBlocks;
  GCBlock *GCBlocks;
} GHMatrix;

// NRow, NColum          - dimension of BlockHMatrix
// NBlockRow,NblockColumn- numbers of HBlocks in row, column
// NDim, MDim            - list starting positions of blocks in
//                         BlockHMatrix
// HBlocks               - list of Hblocks
// Precond ...           - Precondition flag
//                         rest is data specifying the preconditioning
typedef struct sBlockHMatrix
{
  long NRow,NColumn;
  long NBlockRow,NBlockColumn;
  long *NDim, *MDim;

  HMatrix *HBlocks;

  PrecondType Precond;
  long *PrecondIData;
  double complex *PrecondDData; 
  void *PrecondData; 

  long NumBE;
  long *BE;

  csptr csmat;
  iluptr lu;
} BlockHMatrix;

typedef struct sSparseRow
{
  long NI;
  long *J;

  double complex *AIJ; 
} SparseRow;

