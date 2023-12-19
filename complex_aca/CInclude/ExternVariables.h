/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                13.09.2011                                   *
\*****************************************************************************/
/* 
  Last change: 24.11.11 
*/
/*
  Clusters
*/
extern long NCluMin;
extern double ACAEta;
extern double ACAAlpha;
extern double ACABeta;
extern double ACAGamma;
extern long InfoLevel;
/*
  H-Matrix
*/
extern SymmType Symmetry;
extern long MaxRank;
extern double ACAEps;
extern double Memory,AppMemory;
extern PrecondType Precond;
extern long NMaxILU,LFillILU;
extern long LevelBE;
extern long LevelNB;
/*
  Mesh parameters
*/
extern double HMax,HMin,CB;
/*
  Environment
*/
extern char Prefix[4096],Suffix[4096],User[4096];
/*
  Help variables
*/
extern long Bug, Count1, Count2;
/*
  Counter
*/
extern double MemHD,MemHA,MemDH,MemAD,MemAA,MemDD,MemDA,MemAH; 

