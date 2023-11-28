/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                14.11.2010                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HCopyM(CBlock *A,CBlock *B)
{
  HCheckM(A);
  HCheckM(B);
 
  HFreeM(B);

  (*B).IBlock=(*A).IBlock;
  (*B).JBlock=(*A).JBlock;
  (*B).NBlock=(*A).NBlock;
  (*B).MBlock=(*A).MBlock;

  (*B).Type=(*A).Type;

  if((*A).Type  == Null)
    {
      (*B).A11=NULL;
      (*B).A21=NULL;
      (*B).A12=NULL;
      (*B).A22=NULL;

      (*B).Data=NULL;

      (*B).Rank=0;
      (*B).UData=NULL;
      (*B).VData=NULL;
    }
  else if((*A).Type  == Hierarchical)
    {
      (*B).A11=(CBlock *)malloc(1*sizeof(CBlock));
      assert((*B).A11 != NULL);
      HIniM((*B).A11);

      (*B).A21=(CBlock *)malloc(1*sizeof(CBlock));
      assert((*B).A21 != NULL);
      HIniM((*B).A21);

      (*B).A12=(CBlock *)malloc(1*sizeof(CBlock));
      assert((*B).A12 != NULL);
      HIniM((*B).A12);

      (*B).A22=(CBlock *)malloc(1*sizeof(CBlock));
      assert((*B).A22 != NULL);
      HIniM((*B).A22);

      HCopyM((*A).A11,(*B).A11);
      HCopyM((*A).A21,(*B).A21);
      HCopyM((*A).A12,(*B).A12);
      HCopyM((*A).A22,(*B).A22);

      (*B).Data=NULL;

      (*B).Rank=0;
      (*B).UData=NULL;
      (*B).VData=NULL;
    }
  else if ((*A).Type == Dense)
    {
      (*B).A11=NULL;
      (*B).A21=NULL;
      (*B).A12=NULL;
      (*B).A22=NULL;

      (*B).Data=(double complex *)malloc(((*A).NBlock)*((*A).MBlock)*sizeof(double complex));
      assert((*B).Data != NULL);

      cblas_zcopy(((*A).NBlock)*((*A).MBlock),(*A).Data,1,(*B).Data,1);

      (*B).Rank=0;
      (*B).UData=NULL;
      (*B).VData=NULL;
    }
  else if ((*A).Type == Admissible)
    {
      (*B).A11=NULL;
      (*B).A21=NULL;
      (*B).A12=NULL;
      (*B).A22=NULL;

      (*B).Data=NULL;

      (*B).Rank=(*A).Rank;
      (*B).UData=(double complex *)malloc(((*A).NBlock)*((*A).Rank)*sizeof(double complex));
      assert((*B).UData != NULL);

      cblas_zcopy(((*A).NBlock)*((*A).Rank),(*A).UData,1,(*B).UData,1);

      (*B).VData=(double complex *)malloc(((*A).MBlock)*((*A).Rank)*sizeof(double complex));
      assert((*B).VData != NULL);

      cblas_zcopy(((*A).MBlock)*((*A).Rank),(*A).VData,1,(*B).VData,1);
    }

  HCheckM(A);
  HCheckM(B);
}
