/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                13.09.2011                                   *
\*****************************************************************************/

/*
  Main
*/

void Simulation();

/*
  Directory: ACA
*/

void PartialACA(Pair BlockPair,long IBPos,long JBPos, long *PermuRow,long *PermuColumn, Cluster *ClustersRow, Cluster *ClustersColumn, long NBlock,long MBlock, long *aResult, long *aRank,double complex *U,double complex *V, double complex (*ElementName)(long,long), ACAParameters ParamACA);
double complex MatElement(long NBlock,long MBlock,long Rank,double complex *U,double complex *V,long i,long j);
void PartialGACA(GCBlock Block, long *aResult, long *aRank,double complex *U,double complex *V, double complex (*ElementName)(long,long), ACAParameters ParamACA);
void PostCompression(long N,long M,long *aRank, double complex *aU,double complex *aV, double Eps);
/*
  Directory: BasmodC
*/

long lmax(long a,long b);
long lmin(long a,long b);
double max(double a,double b);
double min(double a,double b);
void dset(long N,double A,double *X,long IX);
void zset(long N,double complex A,double complex *X,long IX);
void lset(long N,long A,long *X,long IX);

/*
  Directory: Clusters
*/

long ClusterTree(long N,double *G, double *X,Cluster **aClusters, long **aPermu,ACAParameters ParamACA);
void DevideCluster(double *G,double *X,long *NClusters,long IClu,Cluster **aClusters, long *Permu,ACAParameters ParamACA);
long AdmissiblePair(long IClu1,long IClu2,Cluster *Clusters1,Cluster *Clusters2,ACAParameters ParamACA);
long ClusterPairs(Cluster *ClustersRow,Cluster *ClustersColumn,Pair **aPairs,long *aNMax,long *aMMax,ACAParameters ParamACA);
long SymmClusterPairs(Cluster *Clusters,Pair **aPairs,long *aNMax,long *aMMax,ACAParameters ParamACA);
void PrintCluster(long IClu,Cluster *Clusters);
void PrintClusterInfo(long NClusters,Cluster *Clusters,long *Permu,long NPairs,Pair *Pairs);
long KappaAdmissiblePair(long IClu1,long IClu2, Cluster *Clusters1,Cluster *Clusters2, ACAParameters ParamACA);
long ClusterPoint(long NPoints,Cluster *Clusters,long *Permu,double *Points,GCBlock **aBlocks,long *aNMax,long *aMMax,ACAParameters ParamACA);
void MTVPrintCluster(long N,double *X,Cluster *Clusters,long *Permu,Pair ClPair);
void PrintGClusterInfo(long N,long NClusters,Cluster *Clusters,long *Permu);
void ParaViewPair(long PairNumber,Pair *aPairs,Cluster *aClusters,long *Permu);

/*
  Directory: HLA
*/

void HCopyM(CBlock *A,CBlock *B);
void HGetLU(CBlock *A,double HLUEps);
void HSumMM(CBlock *A,CBlock *B,double HLUEps);
void HMultMM(CBlock *A,CBlock *B,CBlock *C,double HLUEps);
void HScalM(double complex *Alpha,CBlock *A);
void HSolveLULeft(CBlock *A,CBlock *B,double HLUEps);
void HSolveLURight(CBlock *A,CBlock *B,double HLUEps);
void HAToH(CBlock *A,long N11,long M11);
void HDToH(CBlock *A,long N11,long M11);
void HNToD(CBlock *A);
void HAToD(CBlock *A);
void HHToD(CBlock *A);
void HFreeM(CBlock *A);
void HJoinDenseV(CBlock *A,CBlock *A11,CBlock *A12,double HLUEps);
void HJoinDenseH(CBlock *A,CBlock *A11,CBlock *A21,double HLUEps);
void HSplitDenseH(CBlock *A,long N11,CBlock *A11,CBlock *A21);
void HSplitDenseV(CBlock *A,long M11,CBlock *A11,CBlock *A12);
void HInfoM(CBlock *A,HBlockInfo *AInfo);
void HCheckM(CBlock *A);
void HIniM(CBlock *A);
void HPrintBlocksM(HMatrix *aHMat);
void HDToA(CBlock *A,double HLUEps);
double HMemM(CBlock *A);
double HNormM(CBlock *A);
/*
  Directory: HMatrix
*/

void GenHMatrix(long N,long M,long IBPos,long JBPos, long *PermuRow,long *PermuColumn, Cluster *ClustersRow, Cluster *ClustersColumn, long NPairs,Pair *Pairs,HMatrix *aHMat, double complex (*ElementName)(long,long), ACAParameters *ParamACA);
void GenBlock(Pair BlockPair,long IBPos,long JBPos, long *PermuRow,long *PermuColumn, Cluster *ClustersRow, Cluster *ClustersColumn, long IBlock,long JBlock,long NBlock,long MBlock, double complex *Block, double complex (*ElementName)(long,long));
double complex GenElement(long IEl,long JEl);
void BlockHMatrixMult(BlockHMatrix BlockHMat, double complex *X,double complex *Y);
void HMatrixMult(HMatrix HMat,double complex *X,double complex *Y);
void GetHMatrix(long N,double *X,double *G,HMatrix *HMat, double complex (*ElementName)(long,long), ACAParameters *ParamACA);
void FreeBlockHMatrix(BlockHMatrix *aBlockHMat);
void FreeHMatrix(HMatrix *aHMat);
void GetGHMatrix(long N,double *X,double *G,GHMatrix *GHMat, double complex (*ElementName)(long,long), ACAParameters *aParamACA);
void GenGBlock(GCBlock Block, long IBlock,long JBlock,long NBlock,long MBlock, double complex *EBlock, double complex (*ElementName)(long,long));
void GenGHMatrix(long N,long M, long NBlocks,GCBlock *Blocks,GHMatrix *aGHMat, double complex (*ElementName)(long,long), ACAParameters *aParamACA);
void GHMatrixMult(GHMatrix GHMat,double complex *X,double complex *Y);
void GetBlockHMatrix(long *N,double *X,double *G, BlockHMatrix *BlockHMat, double complex (*ElementName)(long,long), ACAParameters *aParamACA);
void SaveHMatrix(HMatrix HMat,char *HMatrixFileName);
/*
  Directory: Input
*/

long ReadInputFile(void);
/*
  Directory: Linker
*/
void maca_solve_fs(double complex *zbglob, int *ndoftot);
void maca_solve(double complex *zbglob, int *ndoftot);
/*
  Directory: Solver
*/
void BlockHMatrixMultZ(double complex *X,double complex *Y,
                       BlockHMatrix *BlockHMat);
void BlockMHMatrixMultZ(double complex *X,double complex *Y,
                        BlockHMatrix *BlockHMat);
void BlockMHMatrixMultZ_fs(double complex *X,double complex *Y,
                        BlockHMatrix *BlockHMat);
void HMatrixMultZ(double complex *X,double complex *Y,HMatrix *HMat);
void GMResPrint(int Iter,double Error);
void GetSparseStructure(HMatrix *HMat);
void WriteSparseMatrix(char *FileName,SparseCOOMatrix SMat);
void GetBlockSparseStructure(BlockHMatrix *BlockHMat);
void BlockPrePrecond(BlockHMatrix *BlockHMat, double complex (*ElementName)(long,long), double complex *y);
void HLUSolveMM(long NRHS,HMatrix *HMat,double complex *X,
                double complex *Y,double HLUEps);

/* int solve_sym_system_cg_mv(long nEq, void (*mv)(double*,double*,void*), void* data, double* rhs, double* solution, double eps, long maxIter) ; */
/* int solve_system_gmres_mv(long nEq, void (*mat_vec)(double*,double*,void*), void* data, void (*pre_con)(double*,double*,void*), void* pdata, double* rhs, double* solution, double eps, long maxIter, void (*prn_fnc)(int,double)) ; */
int zfgmres1(int nn, void (*matvec)(double complex *,double complex *, void *), void * param, double complex *rhs, double complex *sol, double tol, int im, int *itmax, void (*prn_fnc)(int,double));
void SolveTriDiag(long N,double complex *A,double complex *B,double complex *C,double complex *Y,double complex *F);
void PrePrecond(HMatrix *HMat,double complex (*ElementName)(long,long),double complex *y);
