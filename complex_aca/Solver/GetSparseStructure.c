/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                19.09.2010                                   *
\*****************************************************************************/
/*
  Last change 19.09.2010
*/
/*
  Includes for all functions 
*/

# include "Includes.h"

void GetSparseStructure(HMatrix *HMat)
{
/*
  Local variables
*/
  int N,NumNZ,NBlocks,IB,IBlock,JBlock,NBlock,MBlock,Error,IPos;
  int *NumRow,*PosRow,*BegRow;
  int i,j,ib,jb;
  CBlock CB;

  double complex *Values;
  int *Rows,*Columns;

  N=(*HMat).NRow;
  NBlocks=(*HMat).NBlocks;

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
 
  for (IB=0; IB < NBlocks; IB++)
    {
      CB=(*HMat).CBlocks[IB];
/*
  Small block for ILU detected
*/
      if ((int)CB.Type == 1)
	{
	  NBlock=(int)CB.NBlock;
	  MBlock=(int)CB.MBlock;

	  IBlock=(int)CB.IBlock;

          for(i=IBlock; i < IBlock+NBlock; i++) 
	    NumRow[i]+=MBlock;

          NumNZ+=NBlock*MBlock;
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
  for (IB=0; IB < NBlocks; IB++)
    {
      CB=(*HMat).CBlocks[IB];
/*
  Fill small block for ILU
*/
      if ((int)CB.Type == 1)
	{
          NBlock=(int)CB.NBlock;
          MBlock=(int)CB.MBlock;

	  IBlock=(int)CB.IBlock;
	  JBlock=(int)CB.JBlock;

          for(i=0; i < NBlock; i++) 
	    {
	      ib=IBlock+i;

	      for(j=0; j < MBlock; j++)
		{
		  jb=JBlock+j;
                  IPos=BegRow[ib]+PosRow[ib]+j;
		  Values[IPos]=CB.Data[j*NBlock+i];
		  Rows[IPos]=ib;
		  Columns[IPos]=jb;
		}
	      PosRow[ib]+=MBlock;
	    }

	}
    }

  printf("Sparse part is %8.2f MB (%6.2f %%)\n",16.0*((double) NumNZ)/pow(1024.0,2),100.0*((double) NumNZ)/pow(((double) N),2));


  /* FILE *ResultFile; */

  /* ResultFile=fopen("Data.coo","w"); */

  /* for (i=0; i < NumNZ; i++) */
  /*   fprintf(ResultFile,"%8d %8d %25.15e %25.15e\n",Rows[i],Columns[i],creal(Values[i]),cimag(Values[i])); */

  /* fclose(ResultFile); */

  (*HMat).csmat=(csptr)Malloc(sizeof(zSparMat),"GetSparseStructure : zSparMat");
  Error=zCOOcs(N,NumNZ,Values,Columns,Rows,(*HMat).csmat);
  assert(Error == 0);

  int LFil;
  LFil=1;
  (*HMat).lu = (iluptr)Malloc( sizeof(zILUSpar), "GetSparseStructure : zILUSpar" );
  Error=zilukC(LFil,(*HMat).csmat,(*HMat).lu,stderr);
  if (Error != 0) 
    printf("Error in zilukC = %d\n",Error);
  assert(Error == 0);

  free(NumRow); free(PosRow); free(BegRow);
  free(Values); free(Columns); free(Rows);
}

 
