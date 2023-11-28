//------------------------------------------------------------------
//                                                                  
//  routine name  -  maca_init
//                                                                  
//------------------------------------------------------------------
//                                                                  
//  last revision -  Feb 13
//
//  purpose       -  fill,allocate global data structure, i.e.      
//                 1) allocate memory for hierachical matrix        
//                 2) fill the hierarchy-describing arrays           
//                                                                  
//  in:                                                             
//             Nr - number of blocks ... incase you implement
//                  a HBlock version, so far Nr==1 !
//        Ndoftot - dimension of the global quadratic matrix        
//              X - cluster coordinates                             
//              G - cluster weights                                 
//                                                                  
//  out: in GlobalVariables.h  
//                                                                  
//------------------------------------------------------------------

/* Includes for the main function */
# include "../CInclude/GlobalVariables.h"

#define maca_init maca_init_

void maca_init(int *Nr,int *Ndoftot,int *NRELES,int *Nloc,double *X,double *G)
{
  /* local variables */
  long i,j;
  long Ne;
  long NDoftot;
  long NBlocks;

  /* explicit cast */
  NBlocks = (long) *Nr;
  NDoftot = (long) *Ndoftot;
  Ne = (long) *NRELES;

  /* enable Test prints */
  InfoLevel         = 1;

  /* Set global ACA Parameters */
  PARAMACA.NCluMin  = 64;
  PARAMACA.MaxRank  = 500;
  PARAMACA.ACAEta   = 0.8;
  PARAMACA.ACAAlpha = 0.5;
  PARAMACA.ACABeta  = 0.25;
  PARAMACA.ACAGamma = 0.8;
  PARAMACA.ACAEps   = 1.0e-2;
  PARAMACA.Symmetry = ComplexSymmetric; // for EM-scattering problems, for instance
//  PARAMACA.Symmetry = NonSymmetric; // for matrix-valued fundamental solution generating the matrix
  PARAMACA.Info     = 1;

  printf("maca_init: the ACA parameters are ...\n");
  printf("           NCluMin  = %d\n", PARAMACA.NCluMin);
  printf("           MaxRank  = %d\n", PARAMACA.MaxRank);
  printf("           ACAEta   = %3lf\n", PARAMACA.ACAEta);
  printf("           ACAAlpha = %3lf\n", PARAMACA.ACAAlpha);
  printf("           ACABeta  = %3lf\n", PARAMACA.ACABeta);
  printf("           ACAGamma = %3lf\n", PARAMACA.ACAGamma);
  printf("           ACAEps   = %le\n", PARAMACA.ACAEps);
  printf("           Symmetry = %d\n", PARAMACA.Symmetry);

  /* allocate global C-data structure NLoc and copy Fortran's Loc */
  NLOC[0] = (long *)malloc((Ne+1)*sizeof(long));
  for(i=0; i<Ne+1; i++)
  {
    *(NLOC[0]+i) = (long) *(Nloc+i);
  }

  /* construct global cluster tree */
  if (PARAMACA.Info == 1)
  {
    printf("Construct cluster tree ...\n"); 
    fflush(stdout);
  }
  /* mdlt nodes */
  NCLUSTERS[0]=ClusterTree(Ne,G,X,&(CLUSTERS[0]),
                           &(PERMU[0]),PARAMACA);
  printf("  Done with %ld mdle node clusters\n\n",NCLUSTERS[0]);

  /* Construct cluster pairs */
  if (PARAMACA.Info == 1)
  {
    printf("Construct cluster pairs ...\n"); 
    fflush(stdout);
  }
  /* mdlt pairing */
  NPAIRS[0] = SymmClusterPairs(CLUSTERS[0],&PAIRS[0],
                               NMAX+0,MMAX+0,PARAMACA);
  if (PARAMACA.Info == 1)
  {
    printf("  Done with %ld cluster pairs\n",NPAIRS[0]);
    printf("  Maximal cluster point dimensions NMax = %ld, MMax = %ld\n\n",
           NMAX[0],MMAX[0]);
    fflush(stdout);
  }

  /* generate block structure of global matrix */
  if(NBlocks < 2)
  {
    printf("Not a block matrix ! NBlocks = %ld\n",NBlocks);
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
