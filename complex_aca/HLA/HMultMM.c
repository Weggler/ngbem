/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                27.11.2010                                   *
\*****************************************************************************/
/*
  Last change 19.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HMultMM(CBlock *A,CBlock *B,CBlock *C,double HLUEps)
{
  HCheckM(A);
  HCheckM(B);
  HCheckM(C);
/*
  Local variables
*/
  long NA,MA,NB,MB,N11,M11;
  long i,j;
  double complex C0,C1,*F;
  CBlock *D,*E,*G,*B11,*B21,*A11,*A12;
/*
  Constants
*/
  C0=0.0+0.0*I;
  C1=1.0+0.0*I;
/*
  Help blocks
*/
  D=(CBlock *)malloc(1*sizeof(CBlock));
  assert(D != NULL);
  HIniM(D);

  E=(CBlock *)malloc(1*sizeof(CBlock));
  assert(E != NULL);
  HIniM(E);

  G=(CBlock *)malloc(1*sizeof(CBlock));
  assert(G != NULL);
  HIniM(G);

  B11=(CBlock *)malloc(1*sizeof(CBlock));
  assert(B11 != NULL);
  HIniM(B11);

  B21=(CBlock *)malloc(1*sizeof(CBlock));
  assert(B21 != NULL);
  HIniM(B21);

  A11=(CBlock *)malloc(1*sizeof(CBlock));
  assert(A11 != NULL);
  HIniM(A11);

  A12=(CBlock *)malloc(1*sizeof(CBlock));
  assert(A12 != NULL);
  HIniM(A12);
/*
  Dimension of the matrix A
*/
  NA=(*A).NBlock;
  MA=(*A).MBlock;
/*
  Dimension of the matrix B
*/
  NB=(*B).NBlock;
  MB=(*B).MBlock;

  assert(MA == NB);
/*
  Initialisation of the result matrix
*/
  (*C).IBlock=(*A).IBlock;
  (*C).JBlock=(*B).JBlock;
  (*C).NBlock=NA;
  (*C).MBlock=MB;

  if (((*A).Type == Null) || ((*B).Type == Null))
    {
      (*C).Type=Null;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Hierarchical))
    {
      (*C).Type=Hierarchical;

      (*C).A11=(CBlock *)malloc(1*sizeof(CBlock));
      assert((*C).A11 != NULL);
      HIniM((*C).A11);

      (*C).A21=(CBlock *)malloc(1*sizeof(CBlock));
      assert((*C).A21 != NULL);
      HIniM((*C).A21);

      (*C).A12=(CBlock *)malloc(1*sizeof(CBlock));
      assert((*C).A12 != NULL);
      HIniM((*C).A12);

      (*C).A22=(CBlock *)malloc(1*sizeof(CBlock));
      assert((*C).A22 != NULL);
      HIniM((*C).A22);

      (*C).Data=NULL;

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;
/*
  Block 1,1
*/
      HMultMM((*A).A11,(*B).A11,(*C).A11,HLUEps);
      HMultMM((*A).A12,(*B).A21,D,HLUEps);
      HSumMM(D,(*C).A11,HLUEps);
      HFreeM(D);
/*
  Block 2,1
*/
      HMultMM((*A).A21,(*B).A11,(*C).A21,HLUEps);
      HMultMM((*A).A22,(*B).A21,D,HLUEps);
      HSumMM(D,(*C).A21,HLUEps);
      HFreeM(D);
/*
  Block 1,2
*/
      HMultMM((*A).A11,(*B).A12,(*C).A12,HLUEps);
      HMultMM((*A).A12,(*B).A22,D,HLUEps);
      HSumMM(D,(*C).A12,HLUEps);
      HFreeM(D);
/*
  Block 2,2
*/
      HMultMM((*A).A21,(*B).A12,(*C).A22,HLUEps);
      HMultMM((*A).A22,(*B).A22,D,HLUEps);
      HSumMM(D,(*C).A22,HLUEps);
      HFreeM(D);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Dense))
    {
      N11=(*(*A).A11).MBlock;
      HSplitDenseH(B,N11,B11,B21);
      HMultMM((*A).A11,B11,D,HLUEps);
      HMultMM((*A).A12,B21,E,HLUEps);
      HSumMM(D,E,HLUEps);
      HFreeM(D);

      HMultMM((*A).A21,B11,D,HLUEps);
      HMultMM((*A).A22,B21,G,HLUEps);
      HSumMM(D,G,HLUEps);
      HFreeM(D);
      HFreeM(B11); HFreeM(B21);

      HJoinDenseH(C,E,G,HLUEps);

      HFreeM(E); HFreeM(G);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Admissible))
    {
      (*D).IBlock=(*A).IBlock;
      (*D).JBlock=(*B).JBlock;
      (*D).NBlock=NB;
      (*D).MBlock=(*B).Rank;

      (*D).Type=Dense;

      (*D).A11=NULL;
      (*D).A21=NULL;
      (*D).A12=NULL;
      (*D).A22=NULL;

      (*D).Data=(double complex *)malloc(NB*((*B).Rank)*sizeof(double complex));
      assert((*D).Data != NULL);

      (*D).Rank=0;
      (*D).UData=NULL;
      (*D).VData=NULL;

      cblas_zcopy(NB*((*B).Rank),(*B).UData,1,(*D).Data,1);

      HMultMM(A,D,E,HLUEps);

      if((long)(*E).Type == 0)
	HAToD(E);

      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=(*B).Rank;

      (*C).UData=(double complex *)malloc(NA*((*B).Rank)*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(MB*((*B).Rank)*sizeof(double complex));
      assert((*C).VData != NULL);

      cblas_zcopy(NA*((*B).Rank),(*E).Data,1,(*C).UData,1);
      cblas_zcopy(MB*((*B).Rank),(*B).VData,1,(*C).VData,1);

      HFreeM(D); HFreeM(E);
    }   
  else if (((*A).Type == Dense) && ((*B).Type == Hierarchical))
    {
      M11=(*(*B).A11).NBlock;
      HSplitDenseV(A,M11,A11,A12);

      HMultMM(A11,(*B).A11,D,HLUEps);
      HMultMM(A12,(*B).A21,E,HLUEps);
      HSumMM(D,E,HLUEps);
      HFreeM(D);

      HMultMM(A11,(*B).A12,D,HLUEps);
      HMultMM(A12,(*B).A22,G,HLUEps);
      HSumMM(D,G,HLUEps);
      HFreeM(D);
      HFreeM(A11); HFreeM(A12);
 
      HJoinDenseV(C,E,G,HLUEps);

      HFreeM(E); HFreeM(G);
    }
  else if (((*A).Type == Dense) && ((*B).Type == Dense))
    {
      (*C).Type=Dense;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=(double complex *)malloc(NA*MB*sizeof(double complex));
      assert((*C).Data != NULL);

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;
      
      cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,NA,MB,NB,&C1,(*A).Data,NA,(*B).Data,NB,&C0,(*C).Data,NA);
    }
  else if (((*A).Type == Dense) && ((*B).Type == Admissible))
    {
      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=(*B).Rank;
      (*C).UData=(double complex *)malloc(NA*((*B).Rank)*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(MB*((*B).Rank)*sizeof(double complex));
      assert((*C).VData != NULL);

      cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,NA,(*B).Rank,NB,&C1,(*A).Data,NA,(*B).UData,NB,&C0,(*C).UData,NA);
      cblas_zcopy(MB*((*B).Rank),(*B).VData,1,(*C).VData,1);
    }   
  else if (((*A).Type == Admissible) && ((*B).Type== Hierarchical))
    {
      (*D).IBlock=(*A).IBlock;
      (*D).JBlock=(*B).JBlock;
      (*D).NBlock=(*A).Rank;
      (*D).MBlock=MA;

      (*D).Type=Dense;

      (*D).A11=NULL;
      (*D).A21=NULL;
      (*D).A12=NULL;
      (*D).A22=NULL;

      (*D).Data=(double complex *)malloc(((*A).Rank)*MA*sizeof(double complex));
      assert((*D).Data != NULL);

      for (j=0; j < MA; j++)
	for (i=0; i < (*A).Rank ; i++)
	  (*D).Data[j*(*A).Rank+i]=conj((*A).VData[i*MA+j]);

      (*D).Rank=0;
      (*D).UData=NULL;
      (*D).VData=NULL;

      HMultMM(D,B,E,HLUEps);

      if((long)(*E).Type == 0)
	HAToD(E);

      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;


      (*C).UData=(double complex *)malloc(((*A).Rank)*NA*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(((*A).Rank)*MB*sizeof(double complex));
      assert((*C).VData != NULL);

      (*C).Rank=(*A).Rank;
      cblas_zcopy(NA*((*A).Rank),(*A).UData,1,(*C).UData,1);
      for (j=0; j < MB; j++)
       for (i=0; i < (*A).Rank ; i++)
	 (*C).VData[i*MB+j]=conj((*E).Data[j*(*A).Rank+i]);

      HFreeM(D); HFreeM(E);
    }   
  else if (((*A).Type == Admissible) && ((*B).Type == Dense))
    {
      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=(*A).Rank;
      (*C).UData=(double complex *)malloc(NA*((*C).Rank)*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(MB*((*C).Rank)*sizeof(double complex));
      assert((*C).VData != NULL);

      cblas_zcopy(NA*((*A).Rank),(*A).UData,1,(*C).UData,1);
      cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,MB,(*A).Rank,NB,&C1,(*B).Data,NB,(*A).VData,MA,&C0,(*C).VData,MB);
    }   
  else if (((*A).Type == Admissible) && ((*B).Type == Admissible))
    {
      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=lmin((*A).Rank,(*B).Rank);
      (*C).UData=(double complex *)malloc(NA*((*C).Rank)*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(MB*((*C).Rank)*sizeof(double complex));
      assert((*C).VData != NULL);

      F=(double complex *)malloc(((*A).Rank)*((*B).Rank)*sizeof(double complex));
      assert(F != NULL);

      if((*A).Rank <= (*B).Rank)
	{
          cblas_zcopy(NA*((*A).Rank),(*A).UData,1,(*C).UData,1);
	  cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,(*B).Rank,(*A).Rank,NB,&C1,(*B).UData,NB,(*A).VData,MA,&C0,F,(*B).Rank);
	  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,MB,(*A).Rank,(*B).Rank,&C1,(*B).VData,MB,F,(*B).Rank,&C0,(*C).VData,MB);
	}
      else
	{
	  cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,(*A).Rank,(*B).Rank,NB,&C1,(*A).VData,MA,(*B).UData,NB,&C0,F,(*A).Rank);
	  cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,NA,(*B).Rank,(*A).Rank,&C1,(*A).UData,NA,F,(*A).Rank,&C0,(*C).UData,NA);
          cblas_zcopy(MB*((*B).Rank),(*B).VData,1,(*C).VData,1);
	}

      free(F);
    }

  free(D); free(E); free(G); free(B11); free(B21); free(A11); free(A12); 

  HCheckM(A);
  HCheckM(B);   
  HCheckM(C);   
}
