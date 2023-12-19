/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                18.07.2010                                   *
\*****************************************************************************/
/*
  Last change 18.07.2010
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void FreeBlockHMatrix(BlockHMatrix *aBlockHMat)
{
/*
  Local variables
*/
  long i,ib,ic,ir;
  BlockHMatrix BHMat;

  BHMat=*aBlockHMat;
/*
  Dimensions
*/
  free(BHMat.NDim); free(BHMat.MDim);
/*
  All blocks
*/
  ib=0;

  for (ic=0; ic < BHMat.NBlockColumn; ic++)
    {
      for (ir=0; ir < BHMat.NBlockRow; ir++)
	{

	  for (i=0; i < (BHMat.HBlocks[ib]).NBlocks; i++)
	    if ((BHMat.HBlocks[ib]).CBlocks[i].Type == Dense)
	      free((BHMat.HBlocks[ib]).CBlocks[i].Data);
	    else if ((BHMat.HBlocks[ib]).CBlocks[i].Type == Admissible)
	      {
		free((BHMat.HBlocks[ib]).CBlocks[i].UData);
		free((BHMat.HBlocks[ib]).CBlocks[i].VData);
	      }

	  free((BHMat.HBlocks[ib]).CBlocks);

	  if(ic == ir)
	    free((BHMat.HBlocks[ib]).PermuRow);

	  ib++;
	}
    }
/*
  Blocks
*/
  free(BHMat.HBlocks);
}
