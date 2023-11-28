//------------------------------------------------------------------
//
//  routine name         - PartialMACA
//
//------------------------------------------------------------------
//
//  last revision        - Feb 13
//
//  purpose              - generate (cluster-block) matrix-driven 
//                         H-matrix 
//
//  in:
//    BlockPair          - pair (Clu1+Clu2+Type+Aij) 
//    IBPos,JBPos        - pos of the cluster matrix in H-matrix 
//    PermuRow,Column    - permutation of cluster points wrt to
//                         original ordering
//    ClustrsRow,Column  - row cluster, column cluster
//    IBlock,JBlock      - first pos in cluster permuation vector
//    NBlock,MBlock      - number of cluster points
//    NLocRow,Column     - numbers of dofs per cluster point
//    RankClu            - rank of low rank approximation
//    ElementMatrixName  - subroutine generating element matrix
//                         FORTRAN interface
//    ParamACA           - ACA parameters
//    
//  out:
//    U,VT               - matrices determining the approximation
//    Rank               - rank wrt cluster
//    RankClu            - rank of U,V
//    Result            == True  == approximation completed 
//                      == False == approximation not possible
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

/* zgesvd prototype */
void zgesvd_( char* JOBU, char* JOBVT, int* M, int* N, double complex* A, int* LDA, double* S, double complex* Uloc, int* LDU, double complex* VTloc, int* LDVT, double complex* WORK, int* LWORK, double* RWORK, int* INFO);

void PartialMACA(Pair BlockPair,
                 long IBPos,          long JBPos,
                 long *PermuRow,      long *PermuColumn,
                 Cluster *ClustersRow,Cluster *ClustersColumn,
                 long IBlock,         long JBlock,
                 long NBlock,         long MBlock,
                 long *NLocRow,       long *NLocColumn,
                 long *aResult,
                 long *aRankClu,
                 long *aRank,double complex *U,double complex *VT,
                 void (*ElementMatrixName)(long,long,long,long, 
                                           complex double *),
                 ACAParameters ParamACA)
{
/* 
  Local variables 
*/
  int iprint = 0;                   // turn on test prints

  long N0=0, N1=1;                  // dummies fixed values 
                                    // (satisfy valgrind!)
  const long True=1,False=0;        // logical dummies fixed values 

  long i,j,k;                       // running indices
  long II,JJ;                       // ~
  long ii,jj,kk,ll;                 // ~

  long NDim,MDim;                   // effective dimensions 
  int  NLoc,MLoc;                   // dimension of element matrix for 
                                    //   fixed cluster point
  int  Min,Max;                     // for pseudo inverse
  long NMax,MMax;                   // maximal NLoc,MLoc
  long Rank;                        // effective rank
  long RankClu;                     // rank wrt cluster points
  long IPivot=-1,JPivot=-1;         // pivot in row, column
  long IPCol=-1;                    // provisory candidate for next 
                                    //   pivot in row

  long *GenRows,*GenColumns;        // mark generated rows, columns
  double PMax;                      // overall control parameter
  double PRow=0.,PCol=0.;           // save min max singular value wrt 
                                    //   element matrices 
  double MRow=0.,MCol=0.;           // save max max singular value wrt 
                                    //   element matrices
                                    
  double *CSum;                     // control sum

  double FroNorm2;                  // FrobNorm2 squared Frobenius norm 
                                    //   of MACA matrix,
  double Crit;                      // controls approximation accuracy
  double complex *uu,*vvT;          // auxilary matrices to compute FroNorm2 
  double complex trace;             // dummies to compute the error 
                                    //   (for MACA an array)

  long Continue,Restart;            // variables to control the flow
                                    // of the approximation algorithm
  long Result;                      // Result == True if it worked
  double Error;                     // approximation error, see  

  typedef enum {RowColumn=0,ColumnRow=1} CrossTypeType;
  CrossTypeType CrossType;          // MACA driven by row or column, 
                                    //   default: row

/*
  memory for svd by zgesvd 
*/
  char JOBU  = 'A';
  char JOBVT = 'A';
  double* Sloc1;
  double* Sloc2;
  double* RWORK;
  double complex* A;
  double complex* Uloc1;
  double complex* VTloc1;
  double complex* Uloc2;
  double complex* VTloc2;
  double complex* ZWORK;
  int LWORK, INFO;

  /* Memory allocation (depend on cluster dimensions) */
  GenRows=(long *)malloc(NBlock*sizeof(long));
  assert(GenRows != NULL);

  GenColumns=(long *)malloc(MBlock*sizeof(long));
  assert(GenColumns != NULL);

  CSum=(double *)malloc(NBlock*sizeof(double));
  assert(CSum != NULL);

  /* Initialisation */
  Rank=0; RankClu=0; PMax=1.0;

  lset(NBlock,0,GenRows,1);
  lset(MBlock,0,GenColumns,1);
  dset(NBlock,0.0,CSum,1);

  CrossType=RowColumn; // don't change this

  Result=True; Continue=True; Restart=True;

  Error=1.0; Crit=0.0; FroNorm2=0.0;

  /* MACA dimensions, maximal dimension of element matrices */
  NMax = 0; NDim = 0; 
  for(i=0;i<NBlock;i++)
  { 
    NLoc = NLocRow[PermuRow[IBlock+i]+1]
          -NLocRow[PermuRow[IBlock+i]];
    if(NLoc > NMax) NMax = NLoc;
    NDim += NLoc;
  }
  MMax = 0; MDim = 0;
  for(j=0;j<MBlock;j++)
  { 
    MLoc = NLocColumn[PermuColumn[JBlock+j]+1]
          -NLocColumn[PermuColumn[JBlock+j]];
    if(MLoc > MMax) MMax = MLoc;
    MDim += MLoc;
  }

  /* test print */
  if(iprint == 4)
  {
    printf("PartialMACA  IBlock=%ld JBlock=%ld\n", IBlock, JBlock);
    printf("             NBlock=%ld MBlock=%ld\n", NBlock, MBlock);
    printf("             IBPos =%ld JBPos =%ld\n", IBPos,  JBPos);
    printf("             NDim  =%ld MDim  =%ld\n", NDim,   MDim);
    printf("             NMax  =%ld MMax  =%ld\n", NMax,   MMax);
    fflush(stdout);
  }

  /* allocate memory for svd */
  A     = (double complex*) malloc(NMax*MMax*sizeof(double complex));
  Sloc1 = (double*) malloc(MMax*sizeof(double));
  Uloc1 = (double complex*) malloc(NMax*MMax*sizeof(double complex));
  VTloc1= (double complex*) malloc(NMax*MMax*sizeof(double complex));
  Sloc2 = (double*) malloc(MMax*sizeof(double));
  VTloc2= (double complex*) malloc(NMax*MMax*sizeof(double complex));
  Uloc2 = (double complex*) malloc(NMax*MMax*sizeof(double complex));
  LWORK = 5*MMax*NMax;
  RWORK = (double*) malloc(LWORK*sizeof(double));
  ZWORK = (double complex*) malloc(LWORK*sizeof(double complex));

  while ((Continue == True) || (Restart == True))
  {
    if(CrossType == RowColumn)
    {
      if (((Rank+1)*(NDim+MDim) > NDim*MDim) || (RankClu+2 > ParamACA.MaxRank))
      {
        printf("PartialMACA: stop approximation, criterion NOT reached!\n");
        printf("Rank=%ld NDim=%ld MDim=%ld Error=%4.2le\n", 
               Rank,NDim,MDim,Error);
        fflush(stdout);
        Result=False;
        Continue=False;
      }
      /* Check Error */
      else if (Error <= Crit)
      {
//        printf("PartialMACA: stop approximation, criterion reached!\n");
//        fflush(stdout);
        Continue=False;
        Restart=False;

        /* Check column sum */
        i=0;
        while ((Continue == False) && (i < NBlock))
        {
          if( (CSum[i]<1.0e-15*PMax) && (GenRows[i] == 0))
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
        while ((IPivot < NBlock) && (GenRows[IPivot] == 1)) 
        {
          IPivot++;
        }

        /* All rows are generated */
    	if (IPivot == NBlock)
        {
          printf("PartialMACA: all rows generated! \n");
          Continue=False;
          Restart=False;
    	}
        else
        {
          Continue=True;
          Restart=False;
        }
      }

      /* generate cluster point row dim (NLoc x MDim)^t at pos Rank*MDim */
      if(Continue == True)
      {
        GenMBlock(BlockPair,
                  IBPos,       JBPos,
                  PermuRow,    PermuColumn,
                  ClustersRow, ClustersColumn,
                  IPivot,      N0,      // IPivot,0
                  N1,          MBlock,  // 1,MBlock => V = NLoc x MDim 
                  NLocRow,     NLocColumn,
                  VT+Rank*MDim,
                  ElementMatrixName);
        if(iprint == 2)
          printf("IPivot=%3ld\n",IPivot ); 

        /* 1) mark row IPivot as generated
           2) set NLoc == dimension of cluster point matrix */
        GenRows[IPivot]=1;
        NLoc = (int)( NLocRow[PermuRow[IBlock+IPivot]+1]
                     -NLocRow[PermuRow[IBlock+IPivot]]);

        /* check if zero rows have been generated */
        for(ii=0;ii<NLoc;ii++)
        {
          if (cblas_dznrm2(MDim,VT+(Rank+ii)*MDim,1) <= 1.0e-15*PMax)
          {
            printf("PartialMACA: zero row %ld generated at Rank %ld\n",
                   IPivot,Rank); 
            Continue=False;
            Restart=True;
          }
        }	  
      }

      /* run through the cluster point matrices in column
         1) compute the residuum matrices A 
         2) copy the A's on the corresponding position in VT */
      if(Continue == True)
      {
        /* reset PRow */
        PRow = 0.; MRow = 0.;

        /* compute the row residuum */
        JJ = 0;
        II = 0;
        for (i=0;i<IPivot;i++)
          II += (int) ( NLocRow[PermuRow[IBlock+i]+1]
                       -NLocRow[PermuRow[IBlock+i]]);
        for (j=0;j<MBlock;j++)
        {
          MLoc = (int) ( NLocColumn[PermuColumn[JBlock+j]+1]
                        -NLocColumn[PermuColumn[JBlock+j]]);
          for(ii=0;ii<NLoc;ii++)
          {
            for(jj=0;jj<MLoc;jj++) 
            {
              A[ii+jj*MLoc] = conj( VT[(Rank+ii)*MDim+JJ+jj]
                           -MatElement(NDim,MDim,Rank,U,VT,II+ii,JJ+jj));
            }
          }
          for(ii=0;ii<NLoc;ii++)
          {
            for(jj=0;jj<MLoc;jj++) 
            {
              VT[(Rank+jj)*MDim+JJ+ii] = A[ii+jj*MLoc];
            }
          }
          JJ += MLoc;

          /* run through not yet generated columns
             1) compute the svd of the cluster point matrix A, 
             2) save A with maximal minimal singular value Sloc1[NLoc-1]
                and save the corresponding j as JPivot */
          if (GenColumns[j] != 1)
          {
            Min = MLoc;
            Max = NLoc;
            if(NLoc < MLoc) 
            {
              Min = NLoc;
              Max = MLoc;
            }
            dset(MMax,0.0,Sloc1,1);

            zgesvd_(&JOBU,&JOBVT,&NLoc,&MLoc,A,&NLoc,
               Sloc1,Uloc1,&NLoc,VTloc1,&MLoc,ZWORK,&LWORK,RWORK,&INFO);

//            printf("PartialMACA: maximal singular values in Column %ld:
//                    % 3.2lf % 3.2lf\n", j,Sloc1[0],Sloc1[Min-1]); 
//            fflush(stdout);
            if(Sloc1[Min-1] > PRow) 
            {
              JPivot = j;
              PRow = Sloc1[Min-1];
              cblas_dcopy(MLoc,Sloc1,1, Sloc2,1);
              cblas_zcopy(NLoc*NLoc,Uloc1,1, Uloc2,1);
              cblas_zcopy(MLoc*MLoc,VTloc1,1, VTloc2,1);
            }
            if(Sloc1[0] > MRow) 
            {
              MRow = Sloc1[0];
            }
          }
        }
        assert(JJ == MDim);

        /* either stop in case of linear dependence ...
           ... or compute the new block row */
        if (MRow <= 1.0e-15*PMax)
        {
          Continue=False;
        }
        else
        {
          /* update PMax */
          PMax=max(PMax,PRow);

          /* compute A^-1 or A^+ (not tested for pseudo inverse,
             i.e., in case NLoc != MLoc */
          if(NLoc != MLoc)
          {
            printf("PartialMACA:\n");
            printf("Be careful!!!!!!!!!!! Code not tested %d!=%d\n", NLoc,MLoc);
          }
          for(ii=0;ii<NLoc;ii++)
          {
            for(jj=0;jj<MLoc;jj++)
            {
              A[jj+ii*MLoc] = 0.; 
              for(kk=0;kk<Max;kk++) 
              {
                A[jj+ii*MLoc] += conj(VTloc2[kk+ii*MLoc])
                                *conj(Uloc2[jj+kk*NLoc])/Sloc2[kk];
              }
            }
          }

          /* scale residual row VT from right by A^-1 or A^+ */
          JJ = 0.;
          for (j=0;j<MBlock;j++)
          {
            for(ii=0;ii<NLoc;ii++)
            {
              for(jj=0;jj<MLoc;jj++)
              {
                Uloc2[ii+jj*MLoc] = 0.;
                for(kk=0;kk<Max;kk++)
                {
                  Uloc2[ii+jj*MLoc] += VT[(Rank+ii)*MDim+JJ+kk]*A[jj*MLoc+kk];
                }
              }
            }
            for(ii=0;ii<NLoc;ii++)
            {
              for(jj=0;jj<MLoc;jj++)
              {
                VT[(Rank+ii)*MDim+JJ+jj] = Uloc2[jj+ii*MLoc];
                if(iprint == 2) // check on unitary and zero matrices in VT
                {
                  printf("Row %3ld Error=% 6.5lf + I % 6.5lf\n", 
                             (Rank+jj)*MDim+JJ+ii,
                             creal(VT[(Rank+ii)*MDim+JJ+jj]),
                             cimag(VT[(Rank+ii)*MDim+JJ+jj]));
                }
              }
            }
            JJ += MLoc;
          }
          assert(JJ == MDim);
        }
      }

      /* generate matrix column of dim NDim x MLoc at position Rank*NDim */
      if(Continue == True)
      {
        GenMBlock(BlockPair,
                  IBPos,      JBPos,
                  PermuRow,   PermuColumn,
                  ClustersRow,ClustersColumn,
                  N0,         JPivot, // 0,JPivot
                  NBlock,     N1,     // NBlock,1 => NDim x MLoc columnwise
                  NLocRow,    NLocColumn,
                  U+Rank*NDim,
                  ElementMatrixName);
        if(iprint == 2)
          printf("JPivot=%3ld\n",JPivot); 

        /* 1) mark row JPivot as generated
           2) set MLoc == dimension of cluster point matrix */
        GenColumns[JPivot]=1;
        MLoc = (int) ( NLocColumn[PermuColumn[JBlock+JPivot]+1]
                      -NLocColumn[PermuColumn[JBlock+JPivot]]);

        /* check if any zero columns have been generated */
        for(jj=0;jj<MLoc;jj++)
        {
          if (cblas_dznrm2(NDim,U+(Rank+jj)*NDim,1) <= 1.0e-15*PMax)
          {
            printf("PartialMACA: zero col %ld generated at Rank %ld\n",
                   JPivot,Rank); 
          }
        }

        /* reset PCol */
        PCol = 0.; MCol = 0.;

        /* 1) compute the column residuum and update control sum 
           2) update control sum 
           3) find the IPivot condidate in not yet generated rows */
        II = 0;
        JJ = 0;
        for (j=0;j<JPivot;j++)
          JJ += (int) ( NLocColumn[PermuColumn[JBlock+j]+1]
                       -NLocColumn[PermuColumn[JBlock+j]]);
        for (i=0;i<NBlock;i++)
        {
          NLoc = (int)( NLocRow[PermuRow[IBlock+i]+1]
                       -NLocRow[PermuRow[IBlock+i]]);
          for(jj=0; jj<MLoc; jj++)
          {
            for(ii=0; ii<NLoc; ii++)
            {
              A[ii+jj*NLoc] = U[(Rank+jj)*NDim+II+ii]
                             -MatElement(NDim,MDim,Rank,U,VT,II+ii,JJ+jj);
            }
          }
          for(jj=0; jj<MLoc; jj++)
          {
            for(ii=0; ii<NLoc; ii++)
            {
              U[(Rank+jj)*NDim+II+ii] = A[ii+jj*NLoc];

              CSum[i] += cabs(A[ii+jj*NLoc]);

              if(iprint == 2) // check on zero matrices in U
              {
                printf("Col %3ld Error=% 6.5lf + I % 6.5lf\n", 
                           (Rank+jj)*NDim+II+ii,
                           creal(U[(Rank+jj)*NDim+II+ii]),
                           cimag(U[(Rank+jj)*NDim+II+ii]));
              }
            }
          }
          II += NLoc;

          /* find IPivot candidate in not yet generated rows */
          if (GenRows[i] != 1)
          {
            Min = MLoc;
            Max = NLoc;
            if(NLoc < MLoc) 
            {
              Min = NLoc;
              Max = MLoc;
	        }
            dset(MMax,0.0,Sloc1,1);

            /* compute the singular values of matrix svd */
            zgesvd_(&JOBU,&JOBVT,&NLoc,&MLoc,A,&NLoc,
                    Sloc1,Uloc1,&NLoc,VTloc1,&MLoc,ZWORK,&LWORK,RWORK,&INFO);

            /* store the cluster point i with maximal minimal singular value */
//            printf("PartialMACA: maximal singular values in row %ld: 
//                    % 3.2lf % 3.2lf\n", i,Sloc1[0],Sloc1[Min-1]); 
            if(Sloc1[Min-1] > PCol) 
            {
              IPCol = i;
              PCol = Sloc1[Min-1];
            }
            if(Sloc1[0] > MCol) 
            {
              MCol = Sloc1[0];
            }
          } 
        }
        assert(II == NDim);

        /* check linear dependence */
        if (MCol <= 1.0e-15*PMax)
        {
          printf("PartialMACA: zero row candidate %ld generated at Rank %ld\n",
                 IPCol,Rank); 
          Continue=False;
        }
        else
        {
          IPivot=IPCol;
          PMax=max(PMax,PCol);
        }
      }

      /* cross generated, compute the error, 
         note: the Frobenius norm of the Residuum S_{i+1}, i=RankClu-1, is
	  
         S_{i+1} = \sum_{jj=1}^i U_jj V^*_jj + U_{i+1}V^*_{i+1} 
         ||S_{i+1}||_F^2 = 
         = ||S_i||_F^2 
          + tr (\sum_jj V_{i+1} U^*_{i+1} U_jj V^*_jj)
          + tr (\sum_jj V_jj U^*_jj U_{i+1} V^*_{i+1})
          + tr (V_{i+1} U^*_{i+1} U_{i+1} V^*_{i+1})
         = ||S_i||_F^2 +
          + 2 Re (tr\sum_jj (U^*_{i+1} U_jj) (V^*_jj V_{i+1}))
          + tr (U^*_{i+1} U_{i+1}) (V^*_{i+1} V_i)
      */
      if (Continue == True)
      {
        if(NLoc != MLoc)
        {
          printf("PartialMACA:\n");
          printf("Be careful!\n");
          printf("Error calculation does not cover the general case %d!=%d\n",
                 NLoc,MLoc);
        }
        uu  = (double complex*) malloc(MLoc*NLoc*sizeof(double complex));
        vvT = (double complex*) malloc(MLoc*NLoc*sizeof(double complex));
        for (jj=0; jj<=RankClu; jj++) // sum up to RankClu to compute the extra term
        {
          for (j=0; j<NLoc; j++) 
          {
            for (i=0; i<NLoc; i++)
            {
              uu[i+j*NLoc] = 0.0;
              for (k=0; k<NDim; k++) // uu = U_jj U^*_{i+1} 
              {
                uu[i+j*NLoc] += U[(jj*NLoc+i)*NDim+k]*conj(U[(RankClu*NLoc+j)*NDim+k]);
              }
              vvT[i+j*NLoc] = 0.0;
              for (k=0; k<MDim; k++) // vvT = VT^*_jj VT_{i+1} 
              {
                vvT[i+j*NLoc] += conj(VT[(jj*NLoc+i)*MDim+k])*VT[(RankClu*NLoc+j)*MDim+k];
              }
            }
          }
          trace = 0.0;
          for (j=0; j<NLoc; j++) // tr (vvT uu)
          {
            for (k=0; k<NLoc; k++)
            {
               trace += vvT[k+j*NLoc]*uu[k+j*NLoc];
            }
          }
          if(jj < RankClu) 
          {
            FroNorm2 += 2.0*creal(trace);
          }
          else  // the extra term jj==RankClu
          {
            FroNorm2 += creal(trace);
            assert(Error >= 0.0);
            Error     = sqrt(creal(trace));
          }
        }
        if(FroNorm2 < 0.0)		
        {
          printf("\n\nFroNorm2 =%24.16le, trace =%24.16le\n\n", FroNorm2,trace);
          fflush(stdout);
        }
        assert(FroNorm2 >= 0.0);
        Crit = ParamACA.ACAEps*sqrt(FroNorm2);

        if(iprint == 1)
          printf("Crit =%4.2le FroNorm =%24.16le Error =%4.2le\n",
                 Crit,sqrt(FroNorm2),Error);

        /* Update Rank */
        RankClu += 1;
        Rank    += NLoc;
      }
    }
  }
  RankClu += 1; 
  Rank    += NLoc; 

  if(iprint == 1) 
  {
    printf("PartialMACA  admissible block generated ! \n");
    printf("Rank=%ld NDim=%ld MDim=%ld Error=%4.2le\n", 
           Rank,NDim,MDim,Error);
    fflush(stdout);
  }

  /* Free memory */
  *aResult  =Result;
  *aRank    =Rank;
  *aRankClu =RankClu;

  free(A);
  free(Sloc1); free(Uloc1); free(VTloc1); 
  free(Sloc2); free(Uloc2); free(VTloc2); 
  free(RWORK); free(ZWORK);
  free(GenRows); free(GenColumns); free(CSum);

}

//// useful test prints: check explicity on vanishing cross... 
//
//  if(iprint == 3)
//  {
//    double complex app, exa;
//    double complex *test;
//    test = (double complex*) malloc(MDim*NDim*sizeof(double complex));
//    GenMBlock(BlockPair,
//              IBPos,       JBPos,
//              PermuRow,    PermuColumn,
//              ClustersRow, ClustersColumn,
//              N0,          N0,      // 0,0
//              NBlock,      MBlock,  // NBlock,MBlock => NDim x MDim 
//              NLocRow,     NLocColumn,
//              test,
//              ElementMatrixName);
//  
//    for(II=0;II<Rank;II++)
//    {
//      for(JJ=0;JJ<MDim;JJ++)
//      {
//        app = MatElement(NDim,MDim,Rank,U,VT,JJ,II);
//        exa = test[JJ+II*NDim];
//        printf("%2ld: \t % 5.4lf, % 5.4lf  % 5.4lf, % 5.4lf % 5.4lf, % 5.4lf\n", 
//               JJ+II*NDim, creal(exa-app), cimag(exa-app), 
//                           creal(exa),     cimag(exa), 
//                           creal(app),     cimag(app));  
//      }
//    }
//    for(II=0;II<Rank;II++)
//    {
//      for(JJ=0;JJ<MDim;JJ++)
//      {
//        app = MatElement(NDim,MDim,Rank,U,VT,II,JJ);
//        exa = test[II+JJ*NDim];
//        printf("%2ld: \t % 5.4lf, % 5.4lf  % 5.4lf, % 5.4lf % 5.4lf, % 5.4lf\n", 
//               II+JJ*NDim, creal(exa-app), cimag(exa-app), 
//                           creal(exa),     cimag(exa), 
//                           creal(app),     cimag(app));  
//      }
//    }
//    free(test);
//  }

