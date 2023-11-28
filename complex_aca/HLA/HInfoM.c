/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                03.01.2011                                   *
\*****************************************************************************/
/*
  Last change 03.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HInfoM(CBlock *A,HBlockInfo *AInfo)
{
/*
  Local variables
*/
  HBlockInfo *A11Info,*A21Info,*A12Info,*A22Info;
/*
  General
*/
  (*AInfo).IBlock=(*A).IBlock;
  (*AInfo).JBlock=(*A).JBlock;
  (*AInfo).NBlock=(*A).NBlock;
  (*AInfo).MBlock=(*A).MBlock;

  (*AInfo).Type=(*A).Type;
/*
  Subblocks and memory
*/
  if((*A).Type  == Null)
    {
      (*AInfo).NBlocks=1;
      (*AInfo).HBlocks=0;
      (*AInfo).DBlocks=0;
      (*AInfo).ABlocks=0;

      (*AInfo).NMemory=0.0;
      (*AInfo).HMemory=0.0;
      (*AInfo).DMemory=0.0;
      (*AInfo).AMemory=0.0;
    }
  else if((*A).Type  == Hierarchical)
    {
      A11Info=(HBlockInfo *)malloc(1*sizeof(HBlockInfo));
      assert(A11Info != NULL);
      A21Info=(HBlockInfo *)malloc(1*sizeof(HBlockInfo));
      assert(A21Info != NULL);
      A12Info=(HBlockInfo *)malloc(1*sizeof(HBlockInfo));
      assert(A12Info != NULL);
      A22Info=(HBlockInfo *)malloc(1*sizeof(HBlockInfo));
      assert(A22Info != NULL);

      (*AInfo).NBlocks=0;
      (*AInfo).HBlocks=1;
      (*AInfo).DBlocks=0;
      (*AInfo).ABlocks=0;

      (*AInfo).NMemory=0.0;
      (*AInfo).HMemory=4.0*sizeof(CBlock)/(1024.0*1024*0);
      (*AInfo).DMemory=0.0;
      (*AInfo).AMemory=0.0;

      HInfoM((*A).A11,A11Info);
      HInfoM((*A).A21,A21Info);
      HInfoM((*A).A12,A12Info);
      HInfoM((*A).A22,A22Info);

      (*AInfo).NBlocks+=((*A11Info).NBlocks+(*A21Info).NBlocks+(*A12Info).NBlocks+(*A22Info).NBlocks);
      (*AInfo).HBlocks+=((*A11Info).HBlocks+(*A21Info).HBlocks+(*A12Info).HBlocks+(*A22Info).HBlocks);
      (*AInfo).DBlocks+=((*A11Info).DBlocks+(*A21Info).DBlocks+(*A12Info).DBlocks+(*A22Info).DBlocks);
      (*AInfo).ABlocks+=((*A11Info).ABlocks+(*A21Info).ABlocks+(*A12Info).ABlocks+(*A22Info).ABlocks);

      (*AInfo).NMemory+=((*A11Info).NMemory+(*A21Info).NMemory+(*A12Info).NMemory+(*A22Info).NMemory);
      (*AInfo).HMemory+=((*A11Info).NMemory+(*A21Info).NMemory+(*A12Info).NMemory+(*A22Info).NMemory);
      (*AInfo).DMemory+=((*A11Info).DMemory+(*A21Info).DMemory+(*A12Info).DMemory+(*A22Info).DMemory);
      (*AInfo).AMemory+=((*A11Info).AMemory+(*A21Info).AMemory+(*A12Info).AMemory+(*A22Info).AMemory);
    }
  else if ((*A).Type == Dense)
    {
      (*AInfo).NBlocks=0;
      (*AInfo).HBlocks=0;
      (*AInfo).DBlocks=1;
      (*AInfo).ABlocks=0;

      (*AInfo).NMemory=0.0;
      (*AInfo).HMemory=0.0;
      (*AInfo).DMemory=((*A).NBlock)*((*A).MBlock)*sizeof(double complex)/(1024.0*1024*0);;
      (*AInfo).AMemory=0.0;
    }
  else if ((*A).Type == Admissible)
    {
      (*AInfo).NBlocks=0;
      (*AInfo).HBlocks=0;
      (*AInfo).DBlocks=0;
      (*AInfo).ABlocks=1;

      (*AInfo).NMemory=0.0;
      (*AInfo).HMemory=0.0;
      (*AInfo).DMemory=0.0;
      (*AInfo).AMemory=((*A).NBlock+(*A).MBlock)*((*A).Rank)*sizeof(double complex)/(1024.0*1024*0);
    }

  (*AInfo).Memory=(*AInfo).NMemory+(*AInfo).HMemory+(*AInfo).DMemory+(*AInfo).AMemory;  
}
