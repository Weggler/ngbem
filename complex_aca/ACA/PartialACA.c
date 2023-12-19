//------------------------------------------------------------------
//
//  routine name         - PartialMACA
//
//------------------------------------------------------------------
//
//  last revision        - Aug 12
//
//  purpose              - generate cluster-block H-matrix 
//
//  in:
//    BlockPair          - pair (Clu1+Clu2+Type+Aij) 
//    IBPos,JBPos        - position of the cluster matrix in H-matrix 
//    PermuRow,Column    - permutation of cluster points wrt to
//                         original ordering
//    ClustrsRow,Column  - row cluster, column cluster
//    NBlock,MBlock      - number of cluster points
//    Rank               - rank of approximation (rank wrt cluster)
//    ElementName        - subroutine generating full matrix entry
//                         FORTRAN interface
//    ParamACA           - ACA parameters
//    
//  out:
//    U,VT               - matrices determining the approximation
//    Rank               - rank of U,V
//    Result            == True  == approximation completed 
//                      == False == approximation not possible
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

void PartialACA(Pair BlockPair,
                long IBPos,long JBPos,
                long *PermuRow,long *PermuColumn,
                Cluster *ClustersRow, Cluster *ClustersColumn,
                long NBlock,long MBlock,
                long *aResult,
                long *aRank,double complex *U,double complex *V,
                double complex (*ElementName)(long,long),
                ACAParameters ParamACA)
{
/*
  Local variables
*/
  long N0=0, N1=1;
  long i,j;

  long IPivot,JPivot,IPCol;
  double PMax;
  double complex PRow,PCol;
  double complex alpha;

  long Rank;
  long Result;
  double complex Sum,C1,C2,CError;
  double FroNorm2,Crit,Error;

  long *GenRows,*GenColumns;
  double *RSum,*CSum;

  typedef enum {RowColumn=0,ColumnRow=1} CrossTypeType;

  CrossTypeType CrossType;

  long Continue,Restart;

  const long True=1,False=0;
/*
  Memory allocation
*/
  GenRows=(long *)malloc(NBlock*sizeof(long));
  assert(GenRows != NULL);

  GenColumns=(long *)malloc(MBlock*sizeof(long));
  assert(GenColumns != NULL);

  CSum=(double *)malloc(NBlock*sizeof(double));
  assert(CSum != NULL);

  RSum=(double *)malloc(MBlock*sizeof(double));
  assert(RSum != NULL);
/*
  Initialisation
*/
  Rank=0; PMax=1.0;

  lset(NBlock,0,GenRows,1);
  lset(MBlock,0,GenColumns,1);

  dset(NBlock,0.0,CSum,1);
  dset(MBlock,0.0,RSum,1);

  CrossType=RowColumn;

  Result=True;
  Continue=True;
  Restart=True;

  Error=1.0;
  Crit=0.0;
  FroNorm2=0.0;

  while ((Continue == True) || (Restart == True))
  {
    if(CrossType == RowColumn)
    {
      if (((Rank+1)*(NBlock+MBlock) > NBlock*MBlock) || (Rank+2 > ParamACA.MaxRank))
      {
        Result=False;
        Continue=False;
      }
      /* Check Error */
      else if (Error <= Crit)
      {
        Continue=False;
        Restart=False;

        /* Check column sum */
        i=0;
        while ((Continue == False) && (i < NBlock))
        {
          if ((CSum[i] < 1.0e-15*PMax) && (GenRows[i] == 0))
          {
            IPivot=i;
            Continue=True;
            Restart=False;
          }
          i++;
        }
      }

      /* If Restart, look for the first, not yet generated row */
      if (Restart == True)
      {
        Restart=False;
        IPivot=0;
        while ((GenRows[IPivot] == 1) && (IPivot < NBlock)) 
        {
          IPivot++;
        }

        /* All rows are generated */
        if (IPivot == NBlock)
        {
          Continue=False;
          Restart=False;
        }
        else
        {
          Continue=True;
          Restart=False;
        }
      }

      /* Generate cluster point row dim 1 x MBlock at pos Rank*MBlock */
      if(Continue == True)
      {
        GenBlock(BlockPair,
                 IBPos,       JBPos,
                 PermuRow,    PermuColumn,
                 ClustersRow, ClustersColumn,
                 IPivot,      N0,      // IPivot,0
                 N1,          MBlock,  // 1,MBlock => 1 x MBlock
                 V+Rank*MBlock,
                 ElementName);

        /* Mark row IPivot as generated */
        GenRows[IPivot]=1;
        /* Check if the row is zero */
        if ((cblas_dznrm2(MBlock,V+Rank*MBlock,1) <= 1.0e-15*PMax))
        {
          Continue=False;
          Restart=True;
        }
      }	  
/* 
  Compute the residuum and the column pivot position 
*/
      if(Continue == True)
      {
        PRow=0.0+0.0*I;
        for (j=0; j < MBlock; j++)
        {
          /* Update row control sum */
          RSum[j]+=cabs(V[Rank*MBlock+j]);
          /* Residuum wrt row */
          V[Rank*MBlock+j] = conj( V[Rank*MBlock+j]
                                  -MatElement(NBlock,MBlock,Rank,U,V,IPivot,j));
          /* Pivot element */
          if (cabs(V[Rank*MBlock+j]) > cabs(PRow))
          {
            PRow=V[Rank*MBlock+j];
            JPivot=j;
          }
        }
        /* Check linear dependent row */
        if (cabs(PRow) <= 1.0e-15*PMax)
        {
          Continue=False;
        }
        else
        {
          /* Update maximal pivot */
          PMax=max(PMax,cabs(PRow));
          /* Scale the row */
          alpha=(1.0+0.0*I)/PRow;
          cblas_zscal(MBlock,&alpha,V+Rank*MBlock,1);
        }
      }
      /* Generate the column of dim NBlock x 1 at position Rank*NBlock */
      if(Continue == True)
      {
        GenBlock(BlockPair,
                 IBPos,       JBPos,
                 PermuRow,    PermuColumn,
                 ClustersRow, ClustersColumn,
                 N0,          JPivot,
                 NBlock,      N1,
                 U+Rank*NBlock,
                 ElementName);

        /* 1) Compute the column residuum
           2) update control sum 
           3) find the IPivot candidate in not yet generated rows */
        PCol=0.0+0.0*I;
        for (i=0; i < NBlock; i++)
        {
          U[Rank*NBlock+i] = U[Rank*NBlock+i]
                            -MatElement(NBlock,MBlock,Rank,U,V,i,JPivot);

          CSum[i]+=cabs(U[Rank*NBlock+i]);

          if (i != IPivot)
          {
            if (cabs(U[Rank*NBlock+i]) > cabs(PCol))
            {
              PCol=U[Rank*NBlock+i];
              IPCol=i;
            }
          } 
        }
        /* either zero entry (lin dependence) ...
           ... or set new pivot row and update PMax */
        if (cabs(PCol) <= 1.0e-15*PMax)
          Continue=False;
        else
        {
          IPivot=IPCol;
          PMax=max(PMax,cabs(PCol));
        } 
      }
/* 
  Cross generated: compute the error, update rank 
*/
      if (Continue == True)
      {
        Sum=0.0+0.0*I;
        for (i=0; i <= Rank; i++)
        {
          cblas_zdotc_sub(NBlock,U+i*NBlock,1,U+Rank*NBlock,1,&C1);
          cblas_zdotc_sub(MBlock,V+Rank*MBlock,1,V+i*MBlock,1,&C2);
          CError=C1*C2;
          Sum+=CError;
        }
 
        FroNorm2+=2.0*creal(Sum)-creal(CError);
        Crit=ParamACA.ACAEps*sqrt(FroNorm2);
        Error=sqrt(creal(CError));

        /* update rank */
        Rank+=1;
      }
    }
  }
/* 
  Free memory 
*/
  *aResult=Result;
  *aRank=Rank;

  free(GenRows); free(GenColumns); free(RSum); free(CSum);
}
