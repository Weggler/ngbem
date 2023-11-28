/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                11.11.2010                                   *
\*****************************************************************************/
/*
  Last change 14.01.2011
*/
/*
  Includes for all functions
*/

# include "../CInclude/Includes.h"

void HSumMM(CBlock *A,CBlock *B,double HLUEps)
{
  HCheckM(A);
  HCheckM(B);
/*
  Local variables
*/
  long NA,MA,NB,MB,NA11,MA11,NA21,MA21,NA12,MA12,NA22,MA22;
  long j,Rank;
  double complex C0,C1;
  CBlock *C;
/*
  Constants
*/
  C0=0.0+0.0*I;
  C1=1.0+0.0*I;
/*
  Help block
*/
  C=(CBlock *)malloc(1*sizeof(CBlock));
  assert(C != NULL);
  HIniM(C);
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

  assert((NA == NB) && (MA == MB));

  if ((*B).Type == Null)
    {
     HCopyM(A,B);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Hierarchical))
    {
      HSumMM((*A).A11,(*B).A11,HLUEps);
      HSumMM((*A).A21,(*B).A21,HLUEps);
      HSumMM((*A).A12,(*B).A12,HLUEps);
      HSumMM((*A).A22,(*B).A22,HLUEps);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Dense))
    {
      HCopyM(A,C);
      HHToD(C);
      HSumMM(C,B,HLUEps);
      HFreeM(C);
    }
  else if (((*A).Type == Hierarchical) && ((*B).Type == Admissible))
    {
      NA11=(*(*A).A11).NBlock;
      MA11=(*(*A).A11).MBlock;
      
      NA21=(*(*A).A21).NBlock;
      MA21=(*(*A).A21).MBlock;
    
      NA12=(*(*A).A12).NBlock;
      MA12=(*(*A).A12).MBlock;
    
      NA22=(*(*A).A22).NBlock;
      MA22=(*(*A).A22).MBlock;
    
      Rank=(*B).Rank;
    
      if(((NA11+MA11)*Rank < NA11*MA11) &&
         ((NA21+MA21)*Rank < NA21*MA21) &&
     ((NA12+MA12)*Rank < NA12*MA12) &&
     ((NA22+MA22)*Rank < NA22*MA22))
    {
      Count1+=1;
      HAToH(B,NA11,MA11);
    }
      else
    {
      Count2+=1;
      HAToD(B);
    }

      HSumMM(A,B,HLUEps);

      if ((*B).Type == Dense)
	HDToA(B,HLUEps);
    }   
  else if (((*A).Type == Dense) && ((*B).Type == Hierarchical))
    {
/*
  Block 1,1
*/
      NA11=(*(*B).A11).NBlock;
      MA11=(*(*B).A11).MBlock;

      (*C).IBlock=(*A).IBlock;
      (*C).JBlock=(*A).JBlock;
      (*C).NBlock=NA11;
      (*C).MBlock=MA11;

      (*C).Type=Dense;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=(double complex *)malloc(NA11*MA11*sizeof(double complex));
      assert((*C).Data != NULL);

      for (j=0; j < MA11; j++)
	cblas_zcopy(NA11,(*A).Data+j*NA,1,(*C).Data+j*NA11,1);

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;

      HSumMM(C,(*B).A11,HLUEps);

      HFreeM(C);
/*
  Block 2,1
*/
      NA21=(*(*B).A21).NBlock;
      MA21=(*(*B).A21).MBlock;

      (*C).IBlock=(*A).IBlock+NA11;
      (*C).JBlock=(*A).JBlock;
      (*C).NBlock=NA21;
      (*C).MBlock=MA21;

      (*C).Type=Dense;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=(double complex *)malloc(NA21*MA21*sizeof(double complex));
      assert((*C).Data != NULL);

      for (j=0; j < MA21; j++)
	cblas_zcopy(NA21,(*A).Data+j*NA+NA11,1,(*C).Data+j*NA21,1);

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;

      HSumMM(C,(*B).A21,HLUEps);

      HFreeM(C);
/*
  Block 1,2
*/
      NA12=(*(*B).A12).NBlock;
      MA12=(*(*B).A12).MBlock;

      (*C).IBlock=(*A).IBlock;
      (*C).JBlock=(*A).JBlock+MA11;
      (*C).NBlock=NA12;
      (*C).MBlock=MA12;

      (*C).Type=Dense;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=(double complex *)malloc(NA12*MA12*sizeof(double complex));
      assert((*C).Data != NULL);

      for (j=0; j < MA12; j++)
	cblas_zcopy(NA12,(*A).Data+(j+MA11)*NA,1,(*C).Data+j*NA12,1);

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;

      HSumMM(C,(*B).A12,HLUEps);

      HFreeM(C);
/*
  Block 2,2
*/
      NA22=(*(*B).A22).NBlock;
      MA22=(*(*B).A22).MBlock;

      (*C).IBlock=(*A).IBlock+NA11;
      (*C).JBlock=(*A).JBlock+MA11;
      (*C).NBlock=NA22;
      (*C).MBlock=MA22;

      (*C).Type=Dense;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=(double complex *)malloc(NA22*MA22*sizeof(double complex));
      assert((*C).Data != NULL);

      for (j=0; j < MA22; j++)
	cblas_zcopy(NA22,(*A).Data+(j+MA11)*NA+NA11,1,(*C).Data+j*NA22,1);

      (*C).Rank=0;
      (*C).UData=NULL;
      (*C).VData=NULL;

      HSumMM(C,(*B).A22,HLUEps);

      HFreeM(C);
    }   
  else if (((*A).Type == Dense) && ((*B).Type == Dense))
    {
      cblas_zaxpy(NA*MA,&C1,(*A).Data,1,(*B).Data,1);
    }
  else if (((*A).Type == Dense) && ((*B).Type == Admissible))
    {
      (*B).Type=Dense;

      (*B).A11=NULL;
      (*B).A21=NULL;
      (*B).A12=NULL;
      (*B).A22=NULL;

      (*B).Data=(double complex *)malloc(NA*MA*sizeof(double complex));
      assert((*B).Data != NULL);

      cblas_zcopy(NA*MA,(*A).Data,1,(*B).Data,1);
      cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,NB,MB,(*B).Rank,&C1,(*B).UData,NB,(*B).VData,MB,&C1,(*B).Data,NB);

      free((*B).UData); free((*B).VData);

      (*B).Rank=0;
      (*B).UData=NULL;
      (*B).VData=NULL;

      HDToA(B,HLUEps);
    }   
  else if (((*A).Type == Admissible) && ((*B).Type== Hierarchical))
    {
/*
  Block 1,1
*/
      NA11=(*(*B).A11).NBlock;
      MA11=(*(*B).A11).MBlock;

      (*C).IBlock=(*B).IBlock;
      (*C).JBlock=(*B).JBlock;
      (*C).NBlock=NA11;
      (*C).MBlock=MA11;

      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=(*A).Rank;
      (*C).UData=(double complex *)malloc(NA11*(*C).Rank*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(MA11*(*C).Rank*sizeof(double complex));
      assert((*C).VData != NULL);

      for (j=0; j < (*C).Rank; j++)
	{
	 cblas_zcopy(NA11,(*A).UData+j*NB,1,(*C).UData+j*NA11,1);
	 cblas_zcopy(MA11,(*A).VData+j*MB,1,(*C).VData+j*MA11,1);
	}

      HSumMM(C,(*B).A11,HLUEps);

      HFreeM(C);
/*
  Block 2,1
*/
      NA21=(*(*B).A21).NBlock;
      MA21=(*(*B).A21).MBlock;

      (*C).IBlock=(*B).IBlock+NA11;
      (*C).JBlock=(*B).JBlock;
      (*C).NBlock=NA21;
      (*C).MBlock=MA21;

      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=(*A).Rank;
      (*C).UData=(double complex *)malloc(NA21*(*C).Rank*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(MA21*(*C).Rank*sizeof(double complex));
      assert((*C).VData != NULL);

      for (j=0; j < (*C).Rank; j++)
	{
	 cblas_zcopy(NA21,(*A).UData+j*NB+NA11,1,(*C).UData+j*NA21,1);
	 cblas_zcopy(MA21,(*A).VData+j*MB,1,(*C).VData+j*MA21,1);
	}

      HSumMM(C,(*B).A21,HLUEps);

      HFreeM(C);
/*
  Block 1,2
*/
      NA12=(*(*B).A12).NBlock;
      MA12=(*(*B).A12).MBlock;

      (*C).IBlock=(*B).IBlock;
      (*C).JBlock=(*B).JBlock+MA11;
      (*C).NBlock=NA12;
      (*C).MBlock=MA12;

      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=(*A).Rank;
      (*C).UData=(double complex *)malloc(NA12*(*C).Rank*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(MA12*(*C).Rank*sizeof(double complex));
      assert((*C).VData != NULL);

      for (j=0; j < (*C).Rank; j++)
	{
	 cblas_zcopy(NA12,(*A).UData+j*NB,1,(*C).UData+j*NA12,1);
	 cblas_zcopy(MA12,(*A).VData+j*MB+MA11,1,(*C).VData+j*MA12,1);
	}

      HSumMM(C,(*B).A12,HLUEps);

      HFreeM(C);
/*
  Block 2,2
*/
      NA22=(*(*B).A22).NBlock;
      MA22=(*(*B).A22).MBlock;

      (*C).IBlock=(*B).IBlock+NA11;
      (*C).JBlock=(*B).JBlock+MA11;
      (*C).NBlock=NA22;
      (*C).MBlock=MA22;

      (*C).Type=Admissible;

      (*C).A11=NULL;
      (*C).A21=NULL;
      (*C).A12=NULL;
      (*C).A22=NULL;

      (*C).Data=NULL;

      (*C).Rank=(*A).Rank;
      (*C).UData=(double complex *)malloc(NA22*(*C).Rank*sizeof(double complex));
      assert((*C).UData != NULL);
      (*C).VData=(double complex *)malloc(MA22*(*C).Rank*sizeof(double complex));
      assert((*C).VData != NULL);

      for (j=0; j < (*C).Rank; j++)
	{
	 cblas_zcopy(NA22,(*A).UData+j*NB+NA11,1,(*C).UData+j*NA22,1);
	 cblas_zcopy(MA22,(*A).VData+j*MB+MA11,1,(*C).VData+j*MA22,1);
	}

      HSumMM(C,(*B).A22,HLUEps);

      HFreeM(C);
    }   
  else if (((*A).Type == Admissible) && ((*B).Type == Dense))
    {
      cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,NA,MA,(*A).Rank,&C1,(*A).UData,NA,(*A).VData,MA,&C1,(*B).Data,NB);
   }
  else if (((*A).Type == Admissible) && ((*B).Type == Admissible))
    {
      (*B).UData=(double complex *)realloc((*B).UData,NB*((*A).Rank+(*B).Rank)*sizeof(double complex));
      assert((*B).UData != NULL);
      (*B).VData=(double complex *)realloc((*B).VData,MB*((*A).Rank+(*B).Rank)*sizeof(double complex));
      assert((*B).VData != NULL);

      cblas_zcopy(NA*(*A).Rank,(*A).UData,1,(*B).UData+NB*(*B).Rank,1);
      cblas_zcopy(MA*(*A).Rank,(*A).VData,1,(*B).VData+MB*(*B).Rank,1);

      (*B).Rank+=(*A).Rank;

      if ((*B).Rank > lmin(NB,MB))
	{
	  HAToD(B);
	  HDToA(B,HLUEps);
	}
      else
	PostCompression(NB,MB,&((*B).Rank),(*B).UData,(*B).VData,HLUEps);
    }

  free(C);

  HCheckM(A);
  HCheckM(B);
}
