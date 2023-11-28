//------------------------------------------------------------------
//
//  routine name        - GenMHMatrix
//
//------------------------------------------------------------------
//
//  last revision       - Feb 13
//
//  purpose             - generates the matrix-driven HB-matrix  
//
//  in:
//    NLocRow,Column    - list of number of dofs of attributed to 
//                        cluster points in row, column clusters 
//    N,M               - dimension of the HB-matrix
//    IBPos, JBPos      - position of the block in the global HB-Matrix
//    PermuRow,Column   - row, column cluster point permutation 
//    NPairs            - number of pairs to run through ...
//    Pairs             - list of cluster pairs in this block
//    ElementMatrixName - routine generating an element matrix
//                        Note: for MACA, the cluster points are mdle
//                              nodes of variable order !
//    ParamACA          - ACA parameters 
//
//  out: (in GlobalVariable.h)
//    HMat              - HMatrix
//
//------------------------------------------------------------------

/* Includes for all functions */
//# include "../CInclude/Includes.h"
# include "../CInclude/GlobalVariables.h"

void GenMHMatrix(long *NLocRow, long *NLocColumn,
                 long N,long M,
                 long IBPos,long JBPos,
                 long *PermuRow,long *PermuColumn,
                 Cluster *ClustersRow, Cluster *ClustersColumn,
                 long NPairs,Pair *Pairs,
                 HMatrix *aHMat,
                 void (*ElementMatrixName)(long,long,long,long, 
                                           double complex*),
                 ACAParameters *aParamACA)
{
/* 
  Local variables 
*/
  int iprint=0;       // test prints
  long i, j, IPair;   // running indices 
  long N0=0;          // dummies with fixed values (to satisfy valgrind)
  long IBlock,JBlock; // starting of permutation info in PermuRow,Column
  long NBlock,MBlock; // number of points in cluster 
  long NDim,MDim;     // total number of dof for cluster
  long Result,RankClu,Rank; 
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

  /* HMat INTERFACE between cluster data and HMat dimension */
  HMat.NLocRow     = NLocRow; 
  HMat.NLocColumn  = NLocColumn;

  /* HMat cluster data */
  HMat.NPermuRow   = ClustersRow[0].Number; 
  HMat.NPermuColumn= ClustersColumn[0].Number;
  HMat.PermuRow    = PermuRow;
  HMat.PermuColumn = PermuColumn;
  HMat.NBlocks     = NPairs;

  /* HMat properties */
  HMat.Symmetry    = ParamACA.Symmetry;
  HMat.Precond     = Precond;

  /* Test prints */
  if(iprint == 2)
  {
    printf("GenMHMatrix\n");
    printf("HMat.NRow         = %ld\n", HMat.NRow);
    printf("HMat.NColumn      = %ld\n", HMat.NColumn);
    printf("HMat.NPermuRow    = %ld\n", HMat.NPermuRow);
    printf("HMat.NPermuColumn = %ld\n", HMat.NPermuColumn);
    printf("\n");
  }

  /* Memory for the block */
  HMat.CBlocks=(CBlock *)malloc(NPairs*sizeof(CBlock));
  assert(HMat.CBlocks != NULL);
  for (IPair=0; IPair < NPairs; IPair++)
    HIniM(HMat.CBlocks+IPair);

  /* Memory allocation for the low rank block data */
  UBlock=(double complex *)
      malloc(ParamACA.NMax*ParamACA.MaxRank*sizeof(double complex));
//      malloc(ParamACA.NMax*ParamACA.MMax*sizeof(double complex));
  assert(UBlock != NULL);
  VBlock=(double complex *)
      malloc(ParamACA.MMax*ParamACA.MaxRank*sizeof(double complex));
//      malloc(ParamACA.MMax*ParamACA.MMax*sizeof(double complex));
  assert(VBlock != NULL);

  /* Generate blocks */
  ParamACA.Memory   = 0.0;
  ParamACA.AppMemory= 0.0;
  ParamACA.PCMemory = 0.0;
  for (IPair=0; IPair < NPairs; IPair++)
  {
    /* get the type of the paring */
    BlType=Pairs[IPair].Type;

    /* get the position of the permuation in PERMU where the 
       information starts which is relevant for this cluster */
    IBlock=ClustersRow[Pairs[IPair].Clu1].PermuPos;
    JBlock=ClustersColumn[Pairs[IPair].Clu2].PermuPos;

    /* get the number of elements in the clusters */
    NBlock=ClustersRow[Pairs[IPair].Clu1].Number;
    MBlock=ClustersColumn[Pairs[IPair].Clu2].Number;

    /* ...and the resulting dimension depends on order
          of approximation of the mdle nodes
          (cluster point == element == mdle node dof's are
          passed by NLocRow,Column -- pay attention to permutation) */
    NDim = 0; MDim = 0;
    for(i=0;i<NBlock;i++)
    {
      NDim += NLocRow[PermuRow[IBlock+i]+1]
             -NLocRow[PermuRow[IBlock+i]];
    }
    for(i=0;i<MBlock;i++)
    {
      MDim += NLocColumn[PermuColumn[JBlock+i]+1]
             -NLocColumn[PermuColumn[JBlock+i]];
    }

    if(iprint == 2)
    {
      printf("Info for Pairs[%ld]\n", IPair);
      printf("Block Type == %d \n", BlType);
      printf("IBPos  = %2ld \t JBPos  = %2ld\n", IBPos, JBPos);
      printf("IBlock = %2ld \t JBlock = %2ld\n", IBlock, JBlock);
      printf("NBlock = %2ld \t MBlock = %2ld\n", NBlock, MBlock);
      printf("NDim   = %2ld \t MDim   = %2ld\n", NDim, MDim);
    }

    /* copy this info in CBlock data structure */
    HMat.CBlocks[IPair].Type=BlType;

    HMat.CBlocks[IPair].IBlock=IBlock;
    HMat.CBlocks[IPair].JBlock=JBlock;

    HMat.CBlocks[IPair].NBlock=NBlock;
    HMat.CBlocks[IPair].MBlock=MBlock;

    /* Block generation */
    if ((ParamACA.Symmetry == NonSymmetric) || (IBlock <= JBlock))
    {
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

          /* Memory allocation */
          EBlock=(double complex *)
             malloc(NDim*MDim*sizeof(double complex));
          assert(EBlock != NULL);

          /* Generate dense block */
          GenMBlock(Pairs[IPair],
                    IBPos,JBPos,
                    PermuRow,PermuColumn,
                    ClustersRow,ClustersColumn,
                    N0,N0,
                    NBlock,MBlock,
                    NLocRow,NLocColumn,
                    EBlock,
                    ElementMatrixName);
          /* copy the result in data structure arrays */
          HMat.CBlocks[IPair].Data=EBlock;

          /* estimate the memory requirements */
          ParamACA.Memory   +=(double)NDim*
                              (double)MDim*16.0/(1024.0*1024.0);
          ParamACA.AppMemory+=(double)NDim*
                              (double)MDim*16.0/(1024.0*1024.0);
          ParamACA.PCMemory +=(double)NDim*
                              (double)MDim*16.0/(1024.0*1024.0);

          /* nothing else to do */
          HMat.CBlocks[IPair].A11=NULL;
          HMat.CBlocks[IPair].A21=NULL;
          HMat.CBlocks[IPair].A12=NULL;
          HMat.CBlocks[IPair].A22=NULL;
          
          HMat.CBlocks[IPair].Rank=0;
          HMat.CBlocks[IPair].UData=NULL;
          HMat.CBlocks[IPair].VData=NULL;

          break;
        case Admissible:
//        case Dense:

          /* generate low rank MACA */
          PartialMACA(Pairs[IPair],
                      IBPos,JBPos,
                      PermuRow,PermuColumn,
                      ClustersRow,ClustersColumn,
                      IBlock,JBlock,
                      NBlock,MBlock,
                      NLocRow,NLocColumn, 
                      &Result,
                      &RankClu,
                      &Rank,UBlock,VBlock,
                      ElementMatrixName,ParamACA);             
          if (Result == True)
          {
            if (RankClu == 0) 
              HMat.CBlocks[IPair].Type=Null;
            else
            {
              ParamACA.Memory     += 
                    ((double)NDim *(double)MDim)*16.0/(1024.0*1024.0);
              ParamACA.AppMemory += (double)Rank*
                    ((double)NDim +(double)MDim)*16.0/(1024.0*1024.0);
               
//              PostCompression(NDim,MDim,&Rank,UBlock,VBlock,
//                              ParamACA.ACAEps);

              HMat.CBlocks[IPair].RankClu=RankClu;
              HMat.CBlocks[IPair].Rank=Rank;
              
              ParamACA.PCMemory += (double)Rank*
                    ((double)NDim +(double)MDim)*16.0/(1024.0*1024.0);

              /* allocate memory */
              HMat.CBlocks[IPair].UData=(double complex *)
                      malloc(NDim*Rank*sizeof(double complex));
              assert(HMat.CBlocks[IPair].UData != NULL);
              HMat.CBlocks[IPair].VData=(double complex *)
                      malloc(MDim*Rank*sizeof(double complex));
              assert(HMat.CBlocks[IPair].VData != NULL);

              /* copy the result in data structure arrays */
              cblas_zcopy(NDim*Rank,UBlock,1,
                          HMat.CBlocks[IPair].UData,1);
              cblas_zcopy(MDim*Rank,VBlock,1,
                          HMat.CBlocks[IPair].VData,1);

              /* nothing else left to do */
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

            /* allocate memory for dense matrix */
            EBlock=(double complex *)malloc(
                    NDim*MDim*sizeof(double complex));
            assert(EBlock != NULL);

            /* generate dense matrix */
            GenMBlock(Pairs[IPair],
                      IBPos,JBPos,
                      PermuRow,PermuColumn,
                      ClustersRow,ClustersColumn,
                      N0,N0,
                      NBlock,MBlock,
                      NLocRow,NLocColumn,
                      EBlock,
                      ElementMatrixName);

            /* copy the result in data structure arrays */
            HMat.CBlocks[IPair].Data=EBlock;

            HMat.CBlocks[IPair].A11=NULL;
            HMat.CBlocks[IPair].A21=NULL;
            HMat.CBlocks[IPair].A12=NULL;
            HMat.CBlocks[IPair].A22=NULL;
            
            HMat.CBlocks[IPair].Rank=0;
            HMat.CBlocks[IPair].UData=NULL;
            HMat.CBlocks[IPair].VData=NULL;

            /* estimate the memory requirements */
            ParamACA.Memory    +=
                       (double)NDim*
                       (double)MDim*16.0/(1024.0*1024.0);
            ParamACA.AppMemory +=
                       (double)NDim*
                       (double)MDim*16.0/(1024.0*1024.0);
            ParamACA.PCMemory  +=
                       (double)NDim*
                       (double)MDim*16.0/(1024.0*1024.0);
          }          
          break;
      }
    }
  }
  ParamACA.Memory = (double)N*
                    (double)M*16.0/(1024.0*1024.0);

  /* Free memory */
  free(UBlock); 
  free(VBlock);

  *aHMat=HMat;
  *aParamACA=ParamACA;

}

// test
//          int ii,jj;
//          int test[24];
////          double complex *EBlock_new;
//          EBlock_new=(double complex *) malloc(NDim*MDim*sizeof(double complex));
//		  for(i=0;i<NDim;i++)
//          {
//            for(j=0;j<NDim;j++)
//            {
//              EBlock_new[j+i*NDim] = EBlock[j+i*MDim];
//            }
//          }
//          EBlock_new=(double complex *) malloc(144*sizeof(double complex));
//          for(i=0;i<144;i++)
//             EBlock_new[i] = 0.0+0.0*I;
//          test[0] = 0 ; // 2;  0
//          test[1] = 1 ; // 9;  1
//          test[2] = 2 ; // 1;  2
//
//          test[3] = 3 ; // 3;  3
//          test[4] = 4 ; //10;  4
//          test[5] = 0 ; // 2;  0
//          
//		  test[6] = 5 ; // 4;  5
//          test[7] = 6 ; //11;  6
//          test[8] = 3 ; // 3;  3
//          
//		  test[9] = 2 ; // 1;  2
//          test[10] =7 ; //12;  7
//          test[11] =5 ; // 4;  5
//          
//		  test[12] =1 ; // 9;  1
//          test[13] =8 ; // 6;  8
//          test[14] =9 ; // 5;  9
//          
//		  test[15] =4 ; //10;  4
//          test[16] =10; // 7; 10
//          test[17] =8 ; // 6;  8
//          
//		  test[18] =6 ; //11;  6
//          test[19] =11; // 8; 11
//          test[20] =10; // 7; 10
//          
//		  test[21] =7 ; //12;  7
//          test[22] =9 ; // 5;  9
//          test[23] =11; // 8; 11
//          
//		  for(i=0;i<NDim;i++)
//          {
//            ii = test[i]; 
//            for(j=0;j<NDim;j++)
//            {
//              jj = test[j]; 
//              EBlock_new[jj+ii*12] += EBlock[j+i*MDim];
////              EBlock[j+i*MDim] = 0.0+I*0.0;
////              if((j < 24 )&&(i==0))
////                printf("%ld %ld %ld %lf %lf\n",
////                   j+i*MDim,jj,ii,creal(EBlock[i+j*MDim]),creal(EBlock[j+i*MDim]));
////              if((ii == 0 )&&(jj==0))
////              {
////                printf("%ld %ld %ld %lf %lf\n",
////                   j+i*MDim,jj,ii,creal(EBlock[i+j*MDim]),cimag(EBlock[i+j*MDim]));
////                printf("%lf %lf\n",
////                   creal(EBlock_new[jj+ii*12]),cimag(EBlock_new[jj+ii*12]));
////              }
//            }
//          }
//          for(i=0;i<12;i++)
//          {
//            for(j=0;j<12;j++)
//            {
//               EBlock[j+i*12] = EBlock_new[j+i*12];
////               printf("%ld %ld %lf %lf\n",
////                  j,i,creal(EBlock_new[j+i*12]),cimag(EBlock_new[j+i*12]));
//            }
////            printf("\n");
//          }
//          if(iprint == 1) 
//            for(i=0;i<NDim;i++)
//            {
//              for(j=0;j<MDim;j++)
//                 printf("%ld %ld %lf %lf %lf\n",
//                    j,i,cimag(EBlock[j+i*MDim])-cimag(EBlock[i+j*MDim]),
//                        cimag(EBlock[j+i*MDim]),cimag(EBlock[i+j*MDim]));
////                    j,i,creal(EBlock[j+i*MDim])-creal(EBlock[i+j*MDim]),
////                        creal(EBlock[j+i*MDim]),creal(EBlock[i+j*MDim]));
////                    j,i,creal(EBlock[j+i*MDim]),cimag(EBlock[j+i*MDim]));
//              printf("\n");
//            }
//          free(EBlock_new);
