/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                06.04.2009                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void GenGBlock(GCBlock Block,
              long IBlock,long JBlock,long NBlock,long MBlock,
              double complex *EBlock,
              double complex (*ElementName)(long,long))
{
/*
  Local variables
*/
    long IB,JB,IEl,JEl,L,ii;
/*
  Initialisation
*/
    L=0;
/*
  Generation columnwise
*/
    for (JB=JBlock; JB < JBlock+MBlock; JB++)
      {
	JEl=Block.Cols[JB];

	for (IB=IBlock; IB < IBlock+NBlock; IB++)
	  {
	    IEl=Block.Rows[IB];
  if ((IEl < 0) || (IEl >= NGlobal) || (JEl < 0) || (JEl >= NGlobal))
{
  printf("\n"); 
  printf("In BlockGen IEl,JEl= %ld %ld\n",IEl,JEl); 
  printf("   Block N,M= %ld %ld %ld %ld %ld %ld\n",IBlock,JBlock,Block.NBlock,Block.MBlock,NBlock,MBlock);
  printf("     Rows :\n");
  for (ii=0; ii < Block.NBlock; ii++)
    printf(" %ld",Block.Rows[ii]);
  printf("\n");
  printf("     Columns :\n");
  for (ii=0; ii < Block.MBlock; ii++)
    printf(" %ld",Block.Cols[ii]);
  printf("\n");
  getchar();
}

	    EBlock[L]=ElementName(IEl,JEl);

            L+=1;
	  }      
      }
}
