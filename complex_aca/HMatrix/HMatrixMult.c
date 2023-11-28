//------------------------------------------------------------------
//                                                                  
//  routine name  -  HMatrixMult
//                                                                  
//------------------------------------------------------------------
//
//  last revision -  Feb 13
//
//  purpose       -  matrix vector multiplication of a
//                   MH-matrix with a vector
//
//  in:
//           HMat - MH-Matrix
//              X - vector
//
//  out:          
//              Y == HMat X
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

void HMatrixMult(HMatrix HMat,double complex *X,double complex *Y)
{

  /* local variables --- for MACA the local variables do not 
                         coincide with cluster data anymore !!!! */
  int iprint = 0;         // turn on test prints
  long N,M;               // dimension of MH-matrix
  long IB,IE;             // running indices 
  long II,JJ;             // running indices
  long i,j,l;             // running indices
  long IBlock,JBlock;     // for blocks: starting pos of new block contribution.
                          // be careful: in general 
                          // HMat.CBlocks[i].IBlock != IBlock
                          // HMat.CBlocks[i].JBlock != JBlock 
                          // the numbers on the left hold wrt cluster data 
  long NBlock,MBlock;     // dimension of block
                          // be careful: in general 
                          // HMat.CBlocks[i].NBlock != NBlock
                          // HMat.CBlocks[i].MBlock != MBlock 
                          // the numbers on the left hold wrt cluster data 
  long Rank;              // rank == leading dimension  
                          // be careful: in general Rank != RankClu
  double complex *XCopy;  // copies of entry with permutation
  double complex *YCopy;  // temporary result
  double complex *Work;   // library help
  double complex C0,C1;   // Complex zero, complex one
  BlockType BlType;       // needed here to specify the matrix vector multiplication

  /* initialisations */
  C0=0.0+0.0*I;
  C1=1.0+0.0*I;

  /* dimensions of matrix */
  N=HMat.NRow;
  M=HMat.NColumn;

  /* allocate memory for local copies */
  XCopy=(double complex *)malloc(M*sizeof(double complex));
  assert(XCopy != NULL);

  YCopy=(double complex *)malloc(N*sizeof(double complex));
  assert(YCopy != NULL);

  /* copy X taking into account permutation of cluster points */
  if(HMat.NLocColumn == NULL) // == conventional ACA
  {
    for (j=0; j<M; j++)
      XCopy[j]=X[HMat.PermuColumn[j]];
  }
  else // == MACA
  {
    if(iprint == 2) printf("HMatrixMult \n");
    l = 0;
    for (j=0;j<HMat.NPermuColumn;j++)
    {
      IB = HMat.NLocColumn[HMat.PermuColumn[j]];
      IE = HMat.NLocColumn[HMat.PermuColumn[j]+1];
      for(i=IB; i<IE; i++)
      {
        XCopy[l]=X[i];
        l++;
        if(iprint == 1)
        {
          printf("X: %4ld %4ld %4ld %4ld: % 6.5lf+ I % 6.5lf\n",
                 j, HMat.PermuColumn[j], i,l-1,creal(X[i]),cimag(X[i]));
          fflush(stdout);
        }
      }
    }
    assert(l == M);
  }

  zset(N,0.0+0.0*I,YCopy,1);

  /* run through all blocks */
  for (i=0; i<HMat.NBlocks; i++)
  {
     BlType=HMat.CBlocks[i].Type;
  
     if((HMat.NLocColumn == NULL)&&(HMat.NLocRow == NULL)) // == conventional ACA
     {
       IBlock=HMat.CBlocks[i].IBlock;
       JBlock=HMat.CBlocks[i].JBlock;
  
       NBlock=HMat.CBlocks[i].NBlock;
       MBlock=HMat.CBlocks[i].MBlock;
     }
     else // == MACA
     {
       /* hierarchical matrices must not be considered for position ! */
       if ((BlType == Dense)||(BlType == Admissible)) 
       {
         IBlock = 0;
         JBlock = 0;
         for(j=0;j<HMat.CBlocks[i].IBlock;j++)
           IBlock += HMat.NLocRow[HMat.PermuRow[j]+1]
                    -HMat.NLocRow[HMat.PermuRow[j]];
         for(j=0;j<HMat.CBlocks[i].JBlock;j++)
           JBlock += HMat.NLocColumn[HMat.PermuColumn[j]+1]
                    -HMat.NLocColumn[HMat.PermuColumn[j]];
  
         NBlock = 0;
         MBlock = 0;
         for(j=0; j<HMat.CBlocks[i].NBlock; j++)
           NBlock += HMat.NLocRow[HMat.PermuRow[HMat.CBlocks[i].IBlock+j]+1]
                    -HMat.NLocRow[HMat.PermuRow[HMat.CBlocks[i].IBlock+j]];
         for(j=0;j<HMat.CBlocks[i].MBlock;j++)
           MBlock += HMat.NLocColumn[HMat.PermuColumn[HMat.CBlocks[i].JBlock+j]+1]
                    -HMat.NLocColumn[HMat.PermuColumn[HMat.CBlocks[i].JBlock+j]];
       }
     }
     switch(BlType)
     {
       case Dense:
         if((HMat.Symmetry == NonSymmetric) || (IBlock == JBlock))
         {
           cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,MBlock,
                       &C1,HMat.CBlocks[i].Data,NBlock,
                       XCopy+JBlock,1,&C1,YCopy+IBlock,1);
         }
         else if((HMat.Symmetry == ComplexSymmetric) && (IBlock < JBlock))
         {
           cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,MBlock,
                       &C1,HMat.CBlocks[i].Data,NBlock,
                       XCopy+JBlock,1,&C1,YCopy+IBlock,1);
           cblas_zgemv(CblasColMajor,CblasTrans,NBlock,MBlock,
       		           &C1,HMat.CBlocks[i].Data,NBlock,
       		           XCopy+IBlock,1,&C1,YCopy+JBlock,1);
         }
         break;
       case Admissible:
//     if ((BlType == Dense)||(BlType == Admissible))
//     {
         Rank=HMat.CBlocks[i].Rank;
     
         Work=(double complex *)malloc(Rank*sizeof(double complex));
         assert(Work != NULL);

         if((HMat.Symmetry == NonSymmetric) || (IBlock == JBlock))
         {
           cblas_zgemv(CblasColMajor,CblasConjTrans,MBlock,Rank,
                       &C1,HMat.CBlocks[i].VData,MBlock,
                       XCopy+JBlock,1,&C0,Work,1);
     
           cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,Rank,
                       &C1,HMat.CBlocks[i].UData,NBlock,
                       Work,1,&C1,YCopy+IBlock,1);
         }
         else if((HMat.Symmetry == ComplexSymmetric) && (IBlock < JBlock))
         {
           long l;
           double complex *help;
           help=(double complex *)malloc(Rank*MBlock*sizeof(double complex));
           assert(help != NULL);

           cblas_zgemv(CblasColMajor,CblasConjTrans,MBlock,Rank,
                       &C1,HMat.CBlocks[i].VData,MBlock,
                       XCopy+JBlock,1,&C0,Work,1);
     
           cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,Rank,
                       &C1,HMat.CBlocks[i].UData,NBlock,
                       Work,1,&C1,YCopy+IBlock,1);
           
           cblas_zgemv(CblasColMajor,CblasTrans,NBlock,Rank,
                       &C1,HMat.CBlocks[i].UData,NBlock,
                       XCopy+IBlock,1,&C0,Work,1);

           for(l=0; l<Rank*MBlock; l++)
             help[l] = (__real__ HMat.CBlocks[i].VData[l])
                      -I*(__imag__ HMat.CBlocks[i].VData[l]);
         
           cblas_zgemv(CblasColMajor,CblasNoTrans,MBlock,Rank,
                       &C1,help,MBlock,
                       Work,1,&C1,YCopy+JBlock,1);
           free(help);
         }
         free(Work); 
//      }
         break;
    }
  }
     
  if(HMat.NLocRow == NULL)
  {
    for (i=0; i<N; i++)
      Y[HMat.PermuRow[i]]=YCopy[i];
  }
  else
  {
    /* copy result take into account permutation of cluster points */
    if(iprint == 1) printf("HMatrixMult\n");
    l=0;
    for (j=0; j<HMat.NPermuRow; j++)
    {
      IB = HMat.NLocRow[HMat.PermuRow[j]];
      IE = HMat.NLocRow[HMat.PermuRow[j]+1];
      for(i=IB;i<IE;i++)
      {
        Y[i]=YCopy[l];
        l++;
        if(iprint == 1)
        {
           printf("Y: %4ld %4ld  %4ld %4ld: % 6.5lf+ I % 6.5lf\n",
                  j,HMat.PermuRow[j],i,l-1,creal(Y[i]),cimag(Y[i]));
        }
      }
    }
    assert(l == N);
  }

  /* free memory */
  free(XCopy); 
  free(YCopy); 
}

// TEST: 1) generate the full matrix by multiplication of low rank matrics  => test
//       2) perform full matrix vector multiplication by blas routine zgemv => YCopy
//         double complex *test;
//         test = (double complex*) malloc(MBlock*NBlock*sizeof(double complex));
//         for(II=0;II<MBlock;II++)
//           for(JJ=0;JJ<NBlock;JJ++)
//           {
//             test[JJ+II*NBlock] = MatElement(NBlock,MBlock,Rank,
//                                             HMat.CBlocks[i].UData,
//                                             HMat.CBlocks[i].VData,JJ,II);
////             printf("%ld, %4ld % 6.5lf % 6.5lf\n",
////                    i,JJ+II*NBlock,
////                    creal(test[JJ+II*NBlock]),
////                    cimag(test[JJ+II*NBlock]));
//           }
//          cblas_zgemv(CblasColMajor,CblasNoTrans,NBlock,MBlock,
//                      &C1,test,NBlock,
//                      XCopy+JBlock,1,&C1,YCopy+IBlock,1);
//         free(test);

//// ......start svd
//         double complex *zrs, *zls, *zwork;
//         double *rwork, *sv;
//         double smax, smin, cond;
//         int int_info, lwork, ndoftot;
//         lwork = 10*NBlock;
//         sv = (double*) malloc(NBlock*sizeof(double));
//         zrs = (double complex*) malloc(NBlock*NBlock*sizeof(double complex));
//         zls = (double complex*) malloc(NBlock*NBlock*sizeof(double complex));
//         zwork = (double complex*) malloc(lwork*sizeof(double complex));
//         rwork = (double*) malloc(5*NBlock*sizeof(double));
//
//         ndoftot = (int) NBlock;
//         zgesvd_("All","All",&ndoftot,&ndoftot,HMat.CBlocks[i].Data,&ndoftot,
//             sv,zls,&ndoftot,zrs,&ndoftot,zwork,&lwork,rwork,&int_info);
//
//         smax = sv[0];
//         smin = sv[0];
//         for(j=0;j<NBlock;j++)
//         {
//           if(smax < fabs(sv[j]))
//              smax = fabs(sv[j]);
//           if(smin > fabs(sv[j]))
//              smin = fabs(sv[j]);
//           if(j<10)
//             printf("HMatrixMult: %ld, %24.16le \n", j,sv[j]);
//           if(j>NBlock-10)
//             printf("HMatrixMult: %ld, %24.16le \n", j,sv[j]);
//         }
//         cond = smax/smin;
//         printf("HMatrixMult: ndoftot, cond \t %d \t %lf \n", ndoftot, cond);
//         free(zwork);
//         free(rwork);
//         free(sv);
//         free(zls);
//         free(zrs);
////.......end svd

//         for(II=0;II<10;II++)
//           for(JJ=0;JJ<10;JJ++)
//           {
//             printf("%ld %4ld % 6.5lf % 6.5lf\n",
//                    i, JJ+II*NBlock,
//                    creal(HMat.CBlocks[i].Data[JJ+II*NBlock]),
//                    cimag(HMat.CBlocks[i].Data[JJ+II*NBlock]));
//           }
