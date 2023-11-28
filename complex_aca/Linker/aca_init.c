//------------------------------------------------------------------
//
//  routine name  -  aca_init
//
//------------------------------------------------------------------
//
//  last revision -  Aug 12
//
//  purpose       -  fill,allocate global data structure, i.e.
//                 1) allocate memory for block hierachical matrix
//                 2) fill the hierarchy-describing arrays
//
//  in:
//         Nblocks- Number of blocks in row and column (quadratic)
//         Ndoftot- dimension of the global quadratic matrix        
//        NE,NM,NV- Number of Edges, Elements and Vertices
//        XE,XM,XV- cluster coordinates
//        GE,GM,GV- cluster weights
//
//  out: in GlobalVariables.h
//
//------------------------------------------------------------------

/* Includes for the main function */

# include "../CInclude/GlobalVariables.h"

#define aca_init aca_init_
void aca_init(int *Nblocks, 
              int *Ndoftot,
              int    *NE, int    *NM,  
              double *XE, double *XM, 
              double *GE, double *GM)
{
/* local variables */
  long i,j;
  long ne;
  long nm;
  long NBlocks;
  long NDoftot;

/* explicit cast */
  ne = (long) *NE;
  nm = (long) *NM;
  NBlocks = (long) *Nblocks;
  NDoftot = (long) *Ndoftot;

/* enable Test prints */
  InfoLevel         = 1;

/* Set global ACA Parameters */
  PARAMACA.NCluMin  = 64;
  PARAMACA.ACAEta   = 0.8;
  PARAMACA.ACAAlpha = 0.5;
  PARAMACA.ACABeta  = 0.25;
  PARAMACA.ACAGamma = 0.8;
  PARAMACA.ACAEps   = 1.0e-6;
  PARAMACA.Symmetry = ComplexSymmetric; // == ComplexSymmetric works only 
                                        //    for conventional HMatrices
  PARAMACA.MaxRank  = 500;
  PARAMACA.Info     = 1;

  /* construct global cluster trees */
  if (PARAMACA.Info == 1)
  {
    printf("Construct cluster tree ...\n"); 
    fflush(stdout);
  }
  /* ... medg nodes */
  NCLUSTERS[0]=ClusterTree(ne,GE,XE,&(CLUSTERS[0]),
                           &(PERMU[0]),PARAMACA);
  if( nm > 0 )
  {
  /* ... mdlt nodes */
    NCLUSTERS[1]=ClusterTree(nm,GM,XM,&(CLUSTERS[1]),
                             &(PERMU[1]),PARAMACA);
  }
  printf("  Done with %ld medg node clusters\n",  NCLUSTERS[0]);
  printf("  Done with %ld mdle node clusters\n\n",NCLUSTERS[1]);

  /* Construct cluster pairs */
  if (PARAMACA.Info == 1)
  {
    printf("Construct cluster pairs ...\n"); 
    fflush(stdout);
  }
  /* ... medg-medg pairing */
  NPAIRS[0] = SymmClusterPairs(CLUSTERS[0],&PAIRS[0],
                               NMAX+0,MMAX+0,PARAMACA);
  if( nm > 1 )
  {
  /* ...mdlt-mdlt pairing */
    NPAIRS[1] = SymmClusterPairs(CLUSTERS[1],&PAIRS[1],
                                 NMAX+1,MMAX+1,PARAMACA);
  /* ...medg-mdlt pairing */
    NPAIRS[2] =     ClusterPairs(CLUSTERS[0],CLUSTERS[1],&(PAIRS[2]),
                                 NMAX+2,MMAX+2,PARAMACA);
  /* ...mdlt-medg pairing */
    NPAIRS[3] =     ClusterPairs(CLUSTERS[1],CLUSTERS[0],&(PAIRS[3]),
                                 NMAX+3,MMAX+3,PARAMACA);
  }
  for (i=0;i<4;i++)
  {
    if (PARAMACA.Info == 1)
    {
      printf("  Done with %ld cluster pairs\n",NPAIRS[i]);
      printf("  Maximal dimension NMax = %ld, MMax = %ld\n\n",
             NMAX[i],MMAX[i]);
      fflush(stdout);
    }
  }

  /* generate block structure of global matrix */
  if(NBlocks < 2)
  {
    printf("Not a block matrix ! NBlocks = %ld\n",NBlocks);
  }
  else
  {
    printf("A block matrix is generated with NBlocks = %ld\n",NBlocks);
    // change the Symmetry flag to NonSymmetric if necessary
    if(PARAMACA.Symmetry == ComplexSymmetric)
    {
      printf("The complex symmetric case is not implemented yet.\n");
      printf("=> the symmetry flag is changed to NonSymmetric!\n");
      PARAMACA.Symmetry = NonSymmetric;
    }
  }

  /* Memory to store block dimensions */
  /* Block structure and attributes */
  BHMAT = (BlockHMatrix *) malloc(1*sizeof(BlockHMatrix));
  (*BHMAT).Precond = None;
//  (*BHMAT).Precond = BasisElement;
  (*BHMAT).NumBE = 0;
  LevelBE = 4;
  LevelNB = 3;

  (*BHMAT).NRow    = NDoftot;
  (*BHMAT).NColumn = NDoftot;

  (*BHMAT).NBlockRow    = NBlocks;
  (*BHMAT).NBlockColumn = NBlocks;
  printf("init_aca: NBlocks=%ld\n", NBlocks );

  (*BHMAT).NDim=(long *)malloc((*BHMAT).NBlockColumn*sizeof(long));
  assert((*BHMAT).NDim != NULL);

  (*BHMAT).MDim=(long *)malloc((*BHMAT).NBlockColumn*sizeof(long));
  assert((*BHMAT).MDim != NULL);
 
  if (PARAMACA.Info == 1)
  {
    printf("Generate block H-Matrix ...\n"); 
    fflush(stdout);
  }

  (*BHMAT).HBlocks=(HMatrix *)malloc((*BHMAT).NBlockRow*
                                     (*BHMAT).NBlockColumn*sizeof(HMatrix));
  assert((*BHMAT).HBlocks != NULL);

}

