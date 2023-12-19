/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                13.10.2010                                   *
\*****************************************************************************/
/*
  Last change 20.10.2010
*/
/*
  Includes for all functions 
*/

# include "Includes.h"

void GetBlockSparseStructure(BlockHMatrix *BlockHMat)
{
/*
  Local variables
*/
  int N,M,NB,MB,Block,IPos,JPos,IP,MY;
  int NumNZ,NBlocks,IB,IBlock,JBlock,NBlock,MBlock,Error;
  int *NumRow,*PosRow,*BegRow;
  int i,j,ib,jb;
  int LFil;
  CBlock CB;

  double complex *Values;
  int *Rows,*Columns;

  N=(int) (*BlockHMat).NRow;
  M=(int) (*BlockHMat).NColumn;
  
  NB=(int) (*BlockHMat).NBlockRow;
  MB=(int) (*BlockHMat).NBlockColumn;

  NumRow=(int *)malloc(N*sizeof(int));
  assert(NumRow != NULL);

  PosRow=(int *)malloc(N*sizeof(int));
  assert(PosRow != NULL);

  BegRow=(int *)malloc(N*sizeof(int));
  assert(BegRow != NULL);

  for (i=0; i < N; i++)
    {
      NumRow[i]=0;
      PosRow[i]=0;
    }

  NumNZ=0;

  for (jb=0; jb < MB; jb++)
    for (ib=0; ib < NB; ib++)
      {
	IPos=(int) (*BlockHMat).NDim[ib];
 
	Block=jb*NB+ib;

	NBlocks=((*BlockHMat).HBlocks[Block]).NBlocks;

	for (IB=0; IB < NBlocks; IB++)
	  {
	    CB=((*BlockHMat).HBlocks[Block]).CBlocks[IB];
/*
  Small block for ILU detected
*/
	    if ((int) CB.Type == 1)
	      {
		NBlock=(int) CB.NBlock;
		MBlock=(int) CB.MBlock;

		IBlock=(int) CB.IBlock;

		for(i=IPos+IBlock; i < IPos+IBlock+NBlock; i++) 
		  NumRow[i]+=MBlock;

		NumNZ+=NBlock*MBlock;
	      }
	  }
      }
/*
  Memory allocation
*/
  Values=(double complex *)malloc(NumNZ*sizeof(double complex));
  assert(Values != NULL);     

  Rows=(int *)malloc(NumNZ*sizeof(int));
  assert(Rows != NULL);

  Columns=(int *)malloc(NumNZ*sizeof(int));
  assert(Columns != NULL);
/*
  Begin of the rows
*/
  BegRow[0]=0;

  for (i=1; i < N; i++)
    BegRow[i]=BegRow[i-1]+NumRow[i-1];
/*
  Fill values
*/
  for (jb=0; jb < MB; jb++)
    {
      JPos=(int) (*BlockHMat).MDim[jb];

      for (ib=0; ib < NB; ib++)
	{
	  IPos=(int) (*BlockHMat).NDim[ib];
 
	  Block=jb*NB+ib;

	  NBlocks=((*BlockHMat).HBlocks[Block]).NBlocks;

	  for (IB=0; IB < NBlocks; IB++)
	    {
	      CB=((*BlockHMat).HBlocks[Block]).CBlocks[IB];
/*
  Fill small block for ILU
*/
	      if ((int) CB.Type == 1)
		{
		  NBlock=(int) CB.NBlock;
		  MBlock=(int) CB.MBlock;

		  IBlock=(int) CB.IBlock;
		  JBlock=(int) CB.JBlock;

		  for(i=0; i < NBlock; i++) 
		    {
		      ib=IPos+IBlock+i;

		      for(j=0; j < MBlock; j++)
			{
			  jb=JPos+JBlock+j;

			  IP=BegRow[ib]+PosRow[ib]+j;

			  Values[IP]=CB.Data[j*NBlock+i];
			  Rows[IP]=ib;
			  Columns[IP]=jb;
			}

		      PosRow[ib]+=MBlock;
		    }
		}
	    }
	}
    }
/*
  Convert COO format to CS format
*/
  (*BlockHMat).csmat=(csptr)Malloc(sizeof(zSparMat),"GetBlockSparseStructure : zSparMat");
  assert((*BlockHMat).csmat != NULL);

  Error=zCOOcs(N,NumNZ,Values,Columns,Rows,(*BlockHMat).csmat);
  assert(Error == 0);
/*
  Comput ILU decomposition
*/
  (*BlockHMat).lu=(iluptr)Malloc(sizeof(zILUSpar),"GetSparseStructure : zILUSpar");
  assert((*BlockHMat).lu != NULL);

  LFil=(int) LFillILU;
  Error=zilukC(LFil,(*BlockHMat).csmat,(*BlockHMat).lu,stderr);
  assert(Error == 0);

  free(NumRow); free(PosRow); free(BegRow);
  free(Values); free(Columns); free(Rows);
}

 
