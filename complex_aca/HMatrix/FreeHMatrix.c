/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                18.03.2009                                   *
\*****************************************************************************/
/*
  Last change 10.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void FreeHMatrix(HMatrix *aHMat)
{
/*
  Local variables
*/
  long i,CB,CE,i1,i2;
  HMatrix HMat;

  HMat=*aHMat;

  CB=(long)HMat.CBlocks;
  CE=CB+HMat.NBlocks*sizeof(CBlock);
/*
  All blocks
*/
  for (i=0; i < HMat.NBlocks; i++)
    {
      if (HMat.CBlocks[i].Type == Hierarchical)
	{
	  if(((long)HMat.CBlocks[i].A11 < CB) || ((long)HMat.CBlocks[i].A11 > CE))
	  {
	    HFreeM(HMat.CBlocks[i].A11);
	    free(HMat.CBlocks[i].A11);
	    HMat.CBlocks[i].A11=NULL;
	  }
	  if(((long)HMat.CBlocks[i].A21 < CB) || ((long)HMat.CBlocks[i].A21 > CE))
	  {
	    HFreeM(HMat.CBlocks[i].A21);
	    free(HMat.CBlocks[i].A21);
	    HMat.CBlocks[i].A21=NULL;
	  }
	  if(((long)HMat.CBlocks[i].A12 < CB) || ((long)HMat.CBlocks[i].A12 > CE))
	  {
	    HFreeM(HMat.CBlocks[i].A12);
	    free(HMat.CBlocks[i].A12);
	    HMat.CBlocks[i].A12=NULL;
	  }
	  if(((long)HMat.CBlocks[i].A22 < CB) || ((long)HMat.CBlocks[i].A22 > CE))
	  {
	    HFreeM(HMat.CBlocks[i].A22);
	    free(HMat.CBlocks[i].A22);
	    HMat.CBlocks[i].A22=NULL;
	  }
	}
    else if (HMat.CBlocks[i].Type == Dense)
      free(HMat.CBlocks[i].Data);
    else if (HMat.CBlocks[i].Type == Admissible)
      {
	free(HMat.CBlocks[i].UData);
	free(HMat.CBlocks[i].VData);
      }
    }
  free(HMat.CBlocks);

  if (HMat.PermuRow == HMat.PermuColumn)
    free(HMat.PermuRow);
  else 
    {
      free(HMat.PermuRow);
      free(HMat.PermuColumn);
    }

}
