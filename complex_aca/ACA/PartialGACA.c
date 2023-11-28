/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                06.04.2009                                   *
\*****************************************************************************/
/*
  Last change 06.04.2009
*/

/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void PartialGACA(GCBlock Block,
                 long *aResult,
                 long *aRank,double complex *U,double complex *V,
                 double complex (*ElementName)(long,long),
                 ACAParameters ParamACA)
{
/*
  Local variables
*/
  long i,j,Rank,IPivot,JPivot,IPCol,Result,NBlock,MBlock;
  double PMax,FroNorm2,Crit,Error;
  double complex PRow,PCol,alpha,Sum,C1,C2,CError;

  long *GenRows,*GenColumns;
  double *RSum,*CSum;

  typedef enum {RowColumn=0,ColumnRow=1} CrossTypeType;

  CrossTypeType CrossType;

  long Continue,Restart;

  long const True=1,False=0;
/*
  Block size
*/
  NBlock=Block.NBlock;
  MBlock=Block.MBlock;
/*
  Memory allocation
*/
  GenRows=(long *)malloc(NBlock*sizeof(long));
  assert(GenRows != NULL);

  GenColumns=(long *)malloc(MBlock*sizeof(long));
  assert(GenColumns != NULL);

  CSum=(double *)malloc(NBlock*sizeof(double));
  assert(CSum != NULL);

  RSum=(double *)malloc(MBlock*sizeof(double));
  assert(RSum != NULL);
/*
  Initialisation
*/
  Rank=0;
  PMax=1.0;
  FroNorm2=0.0;

  lset(NBlock,0,GenRows,1);
  lset(MBlock,0,GenColumns,1);

  dset(NBlock,0.0,CSum,1);
  dset(MBlock,0.0,RSum,1);
/*
  Sart with Row-Column cross
*/
  CrossType=RowColumn;

  Result=True;
  Continue=True;
  Restart=True;

    while (Continue == True)
      {
	if(CrossType == RowColumn)
	  {
/*
  If Restart, look for the first, not yet generated row
*/
            if (Restart == True)
	      {
		IPivot=0;
		while ((GenRows[IPivot] == 1) && (IPivot < NBlock)) IPivot++;
/*
  All rows are generated
*/
		if (IPivot == NBlock)
		  {
/*
  Free memory
*/
		    *aResult=Result;
		    *aRank=Rank;
		    free(GenRows); free(GenColumns); free(RSum); free(CSum);
		    return;
		  }
	      }
/*
  Generate row
*/
	    GenGBlock(Block,IPivot,0,1,MBlock,V+Rank*MBlock,ElementName);
/*
  Mark row IPivot as generated
*/
	    GenRows[IPivot]=1;
/*
  Zero row detected
*/
	    if (cblas_dznrm2(MBlock,V+Rank*MBlock,1) <= 1.0e-15*PMax)
	      Continue=False;
/*
  Compute the residuum and the column pivot position 
*/
	    else 
	      {
		PRow=0.0+0.0*I;
		for (j=0; j < MBlock; j++)
		  {
/*
  Apdate row control sum
*/
		    RSum[j]+=cabs(V[Rank*MBlock+j]);
/*
  Residuum
*/
		    V[Rank*MBlock+j]=conj(V[Rank*MBlock+j]-MatElement(NBlock,MBlock,Rank,U,V,IPivot,j));
/*
  Pivot element
*/
		    if (cabs(V[Rank*MBlock+j]) > cabs(PRow))
		      {
			PRow=V[Rank*MBlock+j];
			JPivot=j;
		      }
		  }
/*
  Check linear dependent row
*/
		if (cabs(PRow) <= 1.0e-15*PMax)
		  Continue=False;
		else
/*
  Update maximal pivot
*/
		  {
		    PMax=max(PMax,cabs(PRow));
/*
  scale the row
*/
		    alpha=(1.0+0.0*I)/PRow;
		    cblas_zscal(MBlock,&alpha,V+Rank*MBlock,1);
/*
  Generate the column
*/
		    GenGBlock(Block,0,JPivot,NBlock,1,U+Rank*NBlock,ElementName);
/*
  Compute the residuum and the row pivot position 
*/
		    PCol=0.0+0.0*I;
		    for (i=0; i < NBlock; i++)
		      {
/*
  Update column control sum
*/
			CSum[i]+=cabs(U[Rank*NBlock+i]);
/*
  Residuum
*/
			U[Rank*NBlock+i]=U[Rank*NBlock+i]-MatElement(NBlock,MBlock,Rank,U,V,i,JPivot);
/*
  Pivot element
*/
			if (i != IPivot)
			  {
			    if (cabs(U[Rank*NBlock+i]) > cabs(PCol))
			      {
				PCol=U[Rank*NBlock+i];
				IPCol=i;
			      }
			  }
		      }
/*
  Check linear dependent column
*/
		    if (cabs(PCol) <= 1.0e-15*PMax)
		      Continue=False;
		    else
/*
  New pivot row
*/
		      {
			IPivot=IPCol;
/*
  Update maximal pivot
*/
			PMax=max(PMax,cabs(PCol));
		      } 
		  }
	      }
/*
  Compute the error
*/
	    if (Continue == True)
	      {
		Sum=0.0+0.0*I;
		for (i=0; i <= Rank; i++)
		  {
		    cblas_zdotc_sub(NBlock,U+i*NBlock,1,U+Rank*NBlock,1,&C1);
		    cblas_zdotc_sub(MBlock,V+Rank*MBlock,1,V+i*MBlock,1,&C2);
		    CError=C1*C2;
		    Sum+=CError;
		  }

		FroNorm2+=2.0*creal(Sum)-creal(CError);
		Crit=ParamACA.ACAEps*sqrt(FroNorm2);
		Error=sqrt(creal(CError));

		if (Error <= Crit)
		  Continue=False;
                else
                  Restart=False;
	      }
	    if (Continue == False)
	      {
/*
  Control sum check
*/
		i=0;
		Restart=False;
		while ((Restart == False) && (i < NBlock))
		  {
		    if ((CSum[i] < 1.0e-15*PMax) && (GenRows[i] == 0))
		      {
			Rank+=1;
			IPivot=i;
			Restart=True;
			Continue=True;
			printf("Restart forced ! %14.8e  %14.8e\n",CSum[i],PMax);
		      }
		    i++;
		  }
	      }
/*
  Check memory result
*/
	    else
	      {
		if (((Rank+1)*(NBlock+MBlock) > NBlock*MBlock) || (Rank+2 > ParamACA.MaxRank))
		  {
		    Result=False;
		    Continue=False;
		  }
		else
		  Rank+=1;
	      }
	  }
      }
/*
  Free memory
*/
    *aResult=Result;
    *aRank=Rank;

    free(GenRows); free(GenColumns); free(RSum); free(CSum);
}
