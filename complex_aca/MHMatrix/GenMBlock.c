//------------------------------------------------------------------
//
//  routine name         - GenMBlock
//
//------------------------------------------------------------------
//
//  last revision        - Jan 13
//
//  purpose              - generate block matrix-driven H-matrix 
//
//  in:
//    BlockPair          - pair (Clu1+Clu2+Type+Aij)
//    IBPos,JBPos        - position of block matrix in global matrix
//    PermuRow,Column    - permutations from cluster building
//    ClustrsRow,Column  - Clusters 
//    IBlock,JBlock      - cluster row/column numbers to be
//                         generated 
//    NBlock,MBlock      - number of cluster points 
//    ElementMatrixName  - subrouting generating the element matrix
//    NLocRow,Column     - list with dimensions of element matrices
//                         per cluster points
//    
//  out:
//    Block              - a matrix of cluster point dimension 
//                         NBlock x MBlock, the dimension of 
//                         Block is in fact IDim(i)xJDim(j),
//                         0<i<NBlock and 0<J<MBlock
//   0<IMax<NBlock       - (IMax,JMax) specify the ElementMatrix 
//   0<JMax<MBlock         with maximal Norm, where 
//                         Norm = sigma_IDim / sum_i(sigma_i^2), 1<i<IDim
//                         and sigma are the singular values of the 
//                         cluster point element matrix
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

void GenMBlock(Pair BlockPair,
               long IBPos,long JBPos,
               long *PermuRow,long *PermuColumn,
               Cluster *ClustersRow, Cluster *ClustersColumn,
               long IBlock,long JBlock,
               long NBlock,long MBlock,
               long *NLocRow,long *NLocColumn,
               double complex *Block, 
               void (*ElementMatrixName)(long,long,long,long, 
                                         double complex*))
{
/* 
  Explanation of local variables 
*/
  int iprint = 0; // turn on test prints
  long i,j;       // loop through involved cluster points
  long II,JJ;     // positions
  long ii, jj;    // loop through cluster point matrix
  long IEl,JEl;   // real number of cluster points
  long NDim,MDim; // dimension of Block
  long NLoc,MLoc; // dimension of cluster point Matrix
  long IBegin,JBegin; // starting position for the relvant permuation
                      // info, i.e., position in cluster point 
                      // permutation vectors 
  double complex *Matrix;
  double complex help;

/* 
  control input
*/
  if(iprint == 2)
  {
    printf("GenMBlock IBPos =%3ld \tJBPos =%3ld\n", IBPos,JBPos);
    printf("          IBlock=%3ld \tJBlock=%3ld\n", IBlock,JBlock);
    printf("          NBlock=%3ld \tMBlock=%3ld\n", NBlock,MBlock);
    fflush(stdout);
  }
/* 
  Initialisation 
*/
  IBegin=ClustersRow[BlockPair.Clu1].PermuPos   +IBlock; //!= IBlock
  JBegin=ClustersColumn[BlockPair.Clu2].PermuPos+JBlock; //!= JBlock
/* 
  Get the real dimension of Block 
*/
  NDim = 0; 
  for(i=0;i<NBlock;i++) 
  {
    NDim += NLocRow[PermuRow[IBegin+i]+1]
           -NLocRow[PermuRow[IBegin+i]];
  }
  MDim = 0;
  for(j=0;j<MBlock;j++)
  {
    MDim += NLocColumn[PermuColumn[JBegin+j]+1]
           -NLocColumn[PermuColumn[JBegin+j]];
  }
/* 
  control input
*/
  if(iprint == 2)
  {
    printf("          IBegin=%3ld \tJBegin   =%3ld\n", IBegin,JBegin);
    printf("       IPermuPos=%3ld \tJPermuPos=%3ld\n", 
           ClustersRow[BlockPair.Clu1].PermuPos,
           ClustersColumn[BlockPair.Clu2].PermuPos);
    printf("          NDim  =%3ld \tMDim     =%3ld\n\n", NDim,MDim);
    fflush(stdout);
  }
/* 
  Dense Matrix: NBlock > 1 && MBlock > 1 Dimension: NDim x MDim => columnwise stored
  Column      : NBlock > 1 && MBlock = 1 Dimension: NDim x MLoc => U columnwise stored
  Row         : NBlock = 1 && MBlock > 1 Dimension: NLoc x MDim => V store V^t columnwise

  Generate the matrix in three steps
  1) Run through cluster points, MBlock, NBlock
  2) Generate cluster point matrices Matrix, MLoc, NLoc 
  3) Put the cluster point matrices columnwise in Block
*/
  if(NBlock == 1) // store V^t columnwise !!!
  {
    IEl  = PermuRow[IBegin];
    NLoc = NLocRow[IEl+1]-NLocRow[IEl];
    JJ = 0;
    for (j=0; j < MBlock; j++)
    {
      JEl = PermuColumn[JBegin+j];
      MLoc= NLocColumn[JEl+1]-NLocColumn[JEl];
      if(iprint == 2)
      {
        printf("GenMBlock NLoc  =%ld MLoc  =%ld\n", NLoc,MLoc);
        printf("          IEl   =%ld JEl   =%ld\n", IEl,JEl);
        fflush(stdout);
      }

      Matrix = (double complex*) malloc(MLoc*NLoc*
                                        sizeof(double complex));
      for(jj=0;jj<MLoc;jj++)
        for(ii=0;ii<NLoc;ii++)
          Matrix[ii+jj*NLoc] = 0.0;
      ElementMatrixName(IEl,JEl,NLoc,MLoc, Matrix);
      for(jj=0;jj<MLoc;jj++)
      {
        for(ii=0;ii<NLoc;ii++)
        {
           help = *(Matrix+ii+jj*NLoc);
           Block[JJ+jj+ii*MDim] = help;
           if(iprint == 2)
             printf("V^t: %ld (%ld,%ld) from %ld (%ld,%ld) %lf,%lf\n",
                    JJ+jj+ii*MDim,jj+JJ,jj,ii+jj*NLoc,ii,jj,
                    creal(Block[JJ+jj+ii*MDim]),
                    creal(Matrix[ii+jj*NLoc]));
        }
      }
      free(Matrix);
      JJ += MLoc;
    }
    assert(JJ == MDim);
  }
  else
  {
    II = 0;
    JJ = 0;
    for (j=0; j<MBlock; j++)
    {
      JEl = PermuColumn[JBegin+j];
      MLoc= NLocColumn[JEl+1]-NLocColumn[JEl];
      II  = 0;
      for (i=0; i < NBlock; i++)
      {
        IEl =PermuRow[IBegin+i];
//        if((JEl == 0) && (IEl == 1))
//        {
//          printf(" IEl   =%ld JEl   =%ld\n", IEl,JEl);
//          printf(" IBegin=%ld JBegin=%ld\n", IBegin,JBegin);
//          printf(" i=%ld j=%ld\n", i,j);
//        }
        NLoc=NLocRow[IEl+1]-NLocRow[IEl];
        if(iprint == 1)
        {
          printf("GenMBlock NLoc  =%ld MLoc  =%ld\n", NLoc,MLoc);
          printf("          IEl   =%ld JEl   =%ld\n", IEl,JEl);
          fflush(stdout);
        }
        Matrix = (double complex*) malloc(MLoc*NLoc*
                                          sizeof(double complex));

        for(jj=0;jj<MLoc;jj++)
          for(ii=0;ii<NLoc;ii++)
            Matrix[ii+jj*NLoc] = 0.0 + I*0.0;
        ElementMatrixName(IEl,JEl,NLoc,MLoc, Matrix);
        for(jj=0;jj<MLoc;jj++)
        {
          for(ii=0;ii<NLoc;ii++)
          {
            help = *(Matrix+ii+jj*NLoc);
            Block[II+ii+(JJ+jj)*NDim] = help;
            if(iprint == 1)
//            if((JEl == 0) && (IEl == 1))
            {
              printf("U: %ld, (%ld,%ld) from %ld, (%ld,%ld) %lf %lf\n",
                     II+ii+(JJ+jj)*NDim ,II+ii,JJ+jj,ii+jj*NLoc,ii,jj,
                     creal(Block[II+ii+(JJ+jj)*NDim]),
                     creal(Matrix[ii+jj*NLoc]));
            }
          }
        }
        free(Matrix);
        II += NLoc;
      }      
      assert(II == NDim);
      JJ += MLoc;
    }
    assert(JJ == MDim);
  }

}
