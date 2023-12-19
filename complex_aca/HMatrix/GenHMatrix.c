//------------------------------------------------------------------
//
//  routine name        -  GenMHMatrix
//
//------------------------------------------------------------------
//
//  last revision       -  Aug 12
//
//  purpose             -  generates the matrix-driven hierarchic matrix  
//
//  in:
//    N,M               - dimensions of the HMatrix
//    IBPos, JBPos      - position in the global HBlockMatrix
//    PermuRow,Column   - permutation vectors, come from clusters
//    NPairs            - number of cluster pairs
//    Pairs             - list of Pairs
//    ElementMatrixName - routine generating element single (double)
//                        layer potential matrix
//    ParamACA          - ACA parameters 
//
//  out: (in GlobalVariable.h)
//    HMat              - HMatrix
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"
//
//void GenMElement(long IEl,long JEl,long IDim,long JDim, 
//                 double complex *Zloc);

void GenHMatrix(long N,long M,long IBPos,long JBPos,
                long *PermuRow,long *PermuColumn,
                Cluster *ClustersRow, Cluster *ClustersColumn,
                long NPairs,Pair *Pairs,HMatrix *aHMat,
                double complex (*ElementName)(long,long),
                ACAParameters *aParamACA)
{
  /* Local variables */
  long i;
  long IPair,IBlock,JBlock,NBlock,MBlock,Result,Rank,IPara;
  double complex *EBlock,*UBlock,*VBlock;
  HMatrix HMat;
  BlockType BlType;
  ACAParameters ParamACA;
  long const True=1,False=0;

  /* Initialisation */
  ParamACA=*aParamACA;

  /* HMat dimensions */
  HMat.NRow        = N; 
  HMat.NColumn     = M;

  /* HMat interfacing cluster data and HMat dimension -- not needed here */
  HMat.NLocRow     = NULL; 
  HMat.NLocColumn  = NULL;

  /* HMat cluster data */
  HMat.PermuRow    = PermuRow;
  HMat.PermuColumn = PermuColumn;
  HMat.NBlocks     = NPairs;

  /* HMat properties */
  HMat.Symmetry    = ParamACA.Symmetry;
  HMat.Precond     = Precond;
  if(HMat.Symmetry == NonSymmetric)
  {
    printf("GenHMatrix: no symmetry, i.e., HMat.Symmetry = %d\n", 
           HMat.Symmetry);
  }
  else
  {
    printf("GenHMatrix: complex symmetry, i.e., HMat.Symmetry = %d\n", 
           HMat.Symmetry);
  }

  /* Memory for the block */
  HMat.CBlocks=(CBlock *)malloc(NPairs*sizeof(CBlock));
  assert(HMat.CBlocks != NULL);
  for (IPair=0; IPair < NPairs; IPair++)
    HIniM(HMat.CBlocks+IPair);

  /* Memory allocation for the low rank block data */
  UBlock=(double complex *)
      malloc(ParamACA.NMax*ParamACA.MaxRank*sizeof(double complex));
  assert(UBlock != NULL);

  VBlock=(double complex *)
      malloc(ParamACA.MMax*ParamACA.MaxRank*sizeof(double complex));
  assert(VBlock != NULL);

  /* Generate blocks */
  ParamACA.Memory   = 0.0;
  ParamACA.AppMemory= 0.0;
  ParamACA.PCMemory = 0.0;

  IPara=0;

  for (IPair=0; IPair < NPairs; IPair++)
  {
    BlType=Pairs[IPair].Type;
  
    IBlock=ClustersRow[Pairs[IPair].Clu1].PermuPos;
    JBlock=ClustersColumn[Pairs[IPair].Clu2].PermuPos;
  
    NBlock=ClustersRow[Pairs[IPair].Clu1].Number;
    MBlock=ClustersColumn[Pairs[IPair].Clu2].Number;
  
    HMat.CBlocks[IPair].IBlock=IBlock;
    HMat.CBlocks[IPair].JBlock=JBlock;
  
    HMat.CBlocks[IPair].NBlock=NBlock;
    HMat.CBlocks[IPair].MBlock=MBlock;
  
    HMat.CBlocks[IPair].Type=BlType;
  
    /* Block generation */
    switch(BlType)
    {
      case Hierarchical:
        HMat.CBlocks[IPair].A11=HMat.CBlocks+Pairs[IPair].A11;
        HMat.CBlocks[IPair].A21=HMat.CBlocks+Pairs[IPair].A21;
        HMat.CBlocks[IPair].A12=HMat.CBlocks+Pairs[IPair].A12;
        HMat.CBlocks[IPair].A22=HMat.CBlocks+Pairs[IPair].A22;
      
        HMat.CBlocks[IPair].Data=NULL;
      
        HMat.CBlocks[IPair].Rank=0;
        HMat.CBlocks[IPair].UData=NULL;
        HMat.CBlocks[IPair].VData=NULL;

        break;
      case Dense:
        if ((ParamACA.Symmetry == NonSymmetric) || (IBlock <= JBlock))
        {
          /* Memory without approximation */
          ParamACA.Memory+=(double)
                         NBlock*(double)MBlock*16.0/(1024.0*1024.0);
          /* Memory allocation for dense block */
              EBlock=(double complex *)
                 malloc(NBlock*MBlock*sizeof(double complex));
              assert(EBlock != NULL);
          /* Generate dense block */
          GenBlock(Pairs[IPair],IBPos,JBPos,
    	             PermuRow,PermuColumn,
    	             ClustersRow,ClustersColumn,
    	             0,0,NBlock,MBlock,
    	             EBlock,
    	             ElementName);

          HMat.CBlocks[IPair].Data=EBlock;
          
          ParamACA.AppMemory+=(double)
                    NBlock*(double)MBlock*16.0/(1024.0*1024.0);
          ParamACA.PCMemory +=(double)
                    NBlock*(double)MBlock*16.0/(1024.0*1024.0);
          HMat.CBlocks[IPair].A11=NULL;
          HMat.CBlocks[IPair].A21=NULL;
          HMat.CBlocks[IPair].A12=NULL;
          HMat.CBlocks[IPair].A22=NULL;
          
          HMat.CBlocks[IPair].Rank=0;
          HMat.CBlocks[IPair].UData=NULL;
          HMat.CBlocks[IPair].VData=NULL;
        }
        break;
      case Admissible:
        if ((ParamACA.Symmetry == NonSymmetric) || (IBlock <= JBlock))
        {
          /* Memory without approximation */
          ParamACA.Memory+=(double)
                        NBlock*(double)MBlock*16.0/(1024.0*1024.0);
          PartialACA(Pairs[IPair],
                     IBPos,JBPos,
                     PermuRow,PermuColumn,
                     ClustersRow,ClustersColumn,
                     NBlock,MBlock,
                     &Result,
                     &Rank,UBlock,VBlock,
                     ElementName,ParamACA);             
      
          if (Result == True)
          {
            if (Rank == 0) 
              HMat.CBlocks[IPair].Type=Null;
            else
            {
              ParamACA.AppMemory+=(double)Rank*((double)NBlock
                                 +(double)MBlock)*16.0/(1024.0*1024.0);
               
              PostCompression(NBlock,MBlock,&Rank,UBlock,VBlock,
                              ParamACA.ACAEps);

              HMat.CBlocks[IPair].Rank=Rank;
              
              ParamACA.PCMemory+=(double)Rank*((double)NBlock
                                +(double)MBlock)*16.0/(1024.0*1024.0);

              /* Allocation of memory */
              HMat.CBlocks[IPair].UData=(double complex *)
                           malloc(NBlock*Rank*sizeof(double complex));
              assert(HMat.CBlocks[IPair].UData != NULL);
              HMat.CBlocks[IPair].VData=(double complex *)
                           malloc(MBlock*Rank*sizeof(double complex));
              assert(HMat.CBlocks[IPair].VData != NULL);

              /* Copy the result */
              cblas_zcopy(NBlock*Rank,UBlock,1,
                          HMat.CBlocks[IPair].UData,1);
              cblas_zcopy(MBlock*Rank,VBlock,1,
                          HMat.CBlocks[IPair].VData,1);
              
              HMat.CBlocks[IPair].A11=NULL;
              HMat.CBlocks[IPair].A21=NULL;
              HMat.CBlocks[IPair].A12=NULL;
              HMat.CBlocks[IPair].A22=NULL;

              HMat.CBlocks[IPair].Data=NULL;
            }
          }
          else
          {
            HMat.CBlocks[IPair].Type=Dense;
            
            EBlock=(double complex *)
                   malloc(NBlock*MBlock*sizeof(double complex));
            assert(EBlock != NULL);
          
            /* Generate dense block */
            GenBlock(Pairs[IPair],IBPos,JBPos,
                     PermuRow,PermuColumn,
                     ClustersRow,ClustersColumn,
                     0,0,NBlock,MBlock,
                     EBlock,
                     ElementName);
          
             HMat.CBlocks[IPair].Data=EBlock;
             
             HMat.CBlocks[IPair].A11=NULL;
             HMat.CBlocks[IPair].A21=NULL;
             HMat.CBlocks[IPair].A12=NULL;
             HMat.CBlocks[IPair].A22=NULL;
             
             HMat.CBlocks[IPair].Rank=0;
             HMat.CBlocks[IPair].UData=NULL;
             HMat.CBlocks[IPair].VData=NULL;
             
             ParamACA.AppMemory+=(double)NBlock*
                                 (double)MBlock*16.0/(1024.0*1024.0);
             ParamACA.PCMemory +=(double)NBlock*
                                 (double)MBlock*16.0/(1024.0*1024.0);
          }
        }
        break;
    }
    if((ParamACA.Info == 1) && (ParamACA.Memory > 0.0))
    {
      printf(" %6.2f %% done  %8.2f MB instead of %8.2f MB (%6.2f %%)\r", ParamACA.Memory /(((double)N)*((double)M)*16.0/(1024.0*1024.0))*100.0, ParamACA.PCMemory,ParamACA.Memory,ParamACA.PCMemory /ParamACA.Memory*100.0);
      fflush(stdout);
    }
  }

  /* Free memory */
  free(UBlock); 
  free(VBlock);

  *aHMat=HMat;
  *aParamACA=ParamACA;
}
