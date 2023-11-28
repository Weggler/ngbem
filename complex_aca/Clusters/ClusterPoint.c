/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                16.03.2009                                   *
\*****************************************************************************/
/*
  Last change 27.04.2009
*/

/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"
/*
  Local prototypes
*/
void BlockCopy(GCBlock *Blocks,GCBlock **aFinBlocks,long NCheck,long NBlocks);
void SortRowsAngle(double *Points,double *Centre,
                   GCBlock *Blocks,GCBlock **aFinBlocks,
                   long NCheck,long *aNBlocks,double Crit,
                   ACAParameters ParamACA);
void GDevideCluster(double *X,long *aNClusters,long IClu,Cluster **aClusters, 
                    long *Permu,double Crit,ACAParameters ParamACA);
void CheckCand(double *Points,Cluster *Clusters,GCBlock *Blocks,
               long NCheck,ACAParameters ParamACA);
/*
  Functions
*/
long ClusterPoint(long NPoints,Cluster *Clusters,long *Permu,double *Points,GCBlock **aBlocks,long *aNMax,long *aMMax,ACAParameters ParamACA)
{
/*
  Local variables
*/
  long i,AllRows,NCand,NCheck,NBlocks,IBlock,NBlock,MBlock,NMax,MMax,AS1,AS2;
  Pair *Pairs;
  GCBlock *Blocks,*FinBlocks;
  double Crit;
/*
  Initial settings
*/
    NCand=1;
    NCheck=0;
    NBlocks=0;

    Blocks=(GCBlock *)malloc(NCand*sizeof(GCBlock));
    assert(Blocks != NULL);

    Blocks[0].NBlock=NPoints;
    Blocks[0].MBlock=NPoints;

    Blocks[0].Cls=0;

    Blocks[0].Rows=(long *)malloc(Blocks[0].NBlock*sizeof(long));
    assert(Blocks[0].Rows != NULL);

    for (i=0; i < Blocks[0].NBlock; i++)
      Blocks[0].Rows[i]=i;
 
    Blocks[0].Cols=Permu;
/*
  Check all possible clusters
*/
    do
      {



	NBlock=Blocks[NCheck].NBlock;
	MBlock=Blocks[NCheck].MBlock;

/*         printf("%ld %ld %ld %ld %ld\r",NBlocks,NCheck,NCand,NBlock,MBlock); */
/*         fflush(stdout); */

         AllRows=NBlock;

	if ((NBlock <= ParamACA.NCluMin) || (MBlock <= ParamACA.NCluMin))
/*
  Small block
*/
	  {
            Blocks[NCheck].Type=Dense;

            NBlocks+=1;
            BlockCopy(Blocks,&FinBlocks,NCheck,NBlocks);
	  }
	else
/*
  Find admissible rows, distance condition
*/
          {
	    CheckCand(Points,Clusters,Blocks,NCheck,ParamACA);

	    NBlock=Blocks[NCheck].NBlock;
/*
  Realloc rows and give the remaining rows to the children
*/
	    if(NBlock < AllRows)
	      {
/*
  Realloc memory for the candidates
*/
		Blocks=(GCBlock *)realloc(Blocks,(NCand+2)*sizeof(GCBlock));
		assert(Blocks != NULL);
/*
  Fill candidates
*/
		AS1=Clusters[Blocks[NCheck].Cls].Son1;
             
		Blocks[NCand].NBlock=AllRows-NBlock;
		Blocks[NCand].MBlock=Clusters[AS1].Number;
		Blocks[NCand].Cls=AS1;
		Blocks[NCand].Rows=(long *)malloc((AllRows-NBlock)*sizeof(long));
		assert(Blocks[NCand].Rows != NULL);

                for (i=0; i < AllRows-NBlock; i++)
		  Blocks[NCand].Rows[i]=Blocks[NCheck].Rows[NBlock+i]; 
		Blocks[NCand].Cols=Permu+Clusters[AS1].PermuPos;
                
		AS2=Clusters[Blocks[NCheck].Cls].Son2;
             
		Blocks[NCand+1].NBlock=AllRows-NBlock;
		Blocks[NCand+1].MBlock=Clusters[AS2].Number;
		Blocks[NCand+1].Cls=AS2;
		Blocks[NCand+1].Rows=(long *)malloc((AllRows-NBlock)*sizeof(long));
		assert(Blocks[NCand+1].Rows != NULL);

                for (i=0; i < AllRows-NBlock; i++)
		  Blocks[NCand+1].Rows[i]=Blocks[NCheck].Rows[NBlock+i]; 

		Blocks[NCand+1].Cols=Permu+Clusters[AS2].PermuPos;

		NCand+=2;

		if (NBlock > 0)
		  {
		    Blocks[NCheck].Rows=(long *)realloc(Blocks[NCheck].Rows,NBlock*sizeof(long));
                    assert(Blocks[NCheck].Rows != NULL);
		  }
	      }
/*
  Small block
*/
	    if (NBlock <= ParamACA.NCluMin)
	      {
		if (NBlock > 0)
		  {
		    Blocks[NCheck].Type=Dense;

		    NBlocks+=1;
		    BlockCopy(Blocks,&FinBlocks,NCheck,NBlocks);
		  }
	      }
	    else
/*
  sort  admissible rows, angle condition
*/
	      {
                Crit=M_PI/(1.0+ParamACA.ACAGamma*Kappa*Clusters[Blocks[NCheck].Cls].Radius);
		SortRowsAngle(Points,Clusters[Blocks[NCheck].Cls].Centre,
                              Blocks,&FinBlocks,NCheck,&NBlocks,Crit,ParamACA);
	      }

          }
	NCheck+=1;
      }
    while (NCheck < NCand);
/*
  Maximal block size
*/
    NMax=1;
    MMax=1;

    for (IBlock=0; IBlock < NBlocks; IBlock++)
      { 
/*         printf("%ld %ld %ld\n",IBlock,FinBlocks[IBlock].NBlock,FinBlocks[IBlock].MBlock); */
/*         getchar(); */
	NMax=lmax(FinBlocks[IBlock].NBlock,NMax);
	MMax=lmax(FinBlocks[IBlock].MBlock,MMax);
      }

    *aBlocks=FinBlocks;
    *aNMax=NMax;
    *aMMax=MMax;

    return NBlocks;    
}
/*
  Block copy
*/
void BlockCopy(GCBlock *Blocks,GCBlock **aFinBlocks,long NCheck,long NBlocks)
{
  long j;

  if (NBlocks == 1)
    {
      (*aFinBlocks)=(GCBlock *)malloc(NBlocks*sizeof(GCBlock));
      assert((*aFinBlocks) != NULL);
    }
  else
    {
      (*aFinBlocks)=(GCBlock *)realloc((*aFinBlocks),NBlocks*sizeof(GCBlock));
      assert((*aFinBlocks) != NULL);
    }
/*
  Copy
*/
  (*aFinBlocks)[NBlocks-1]=Blocks[NCheck];
}
/*
  Sort raws corresponding to the angle
*/
void SortRowsAngle(double *Points,double *Centre,
                   GCBlock *Blocks,GCBlock **aFinBlocks,
                   long NCheck,long *aNBlocks,double Crit,
                   ACAParameters ParamACA)
{
/*
  Local variables
*/
  long N,i,j,IClu,NClusters,*Permu,NBlocks;
  Cluster *Clusters;
  double *X,alpha;

  NBlocks=*aNBlocks;
/*
  Number of points
*/
  N=Blocks[NCheck].NBlock;
/*
  Memory
*/
  X=(double *)malloc(3*N*sizeof(double));
  assert(X != NULL);

  Permu=(long *)malloc(N*sizeof(long));
  assert(Permu != NULL);
/*
  Copy and projection, initial permutation
*/
  for (i=0; i < N; i++)
    {
      cblas_dcopy(3,Points+3*Blocks[NCheck].Rows[i],1,X+3*i,1);
      cblas_daxpy(3,-1.0,Centre,1,X+3*i,1);
      alpha=1.0/cblas_dnrm2(3,X+3*i,1);
      cblas_dscal(3,alpha,X+3*i,1); 
      Permu[i]=i;
   }
/*
  Initial cluster
*/
  NClusters=1;
  IClu=0;

  Clusters=(Cluster *)malloc(NClusters*sizeof(Cluster));
  assert(Clusters != NULL);
/*
  Initial data
*/
  Clusters[IClu].Level=0;
  Clusters[IClu].Father=-1;
  Clusters[IClu].Number=N;
  Clusters[IClu].PermuPos=0;
/*
  Devide cluster
*/
  do 
    {
      GDevideCluster(X,&NClusters,IClu,&Clusters,Permu,Crit,ParamACA);
      IClu+=1;
    }
  while (IClu < NClusters);
/*
  Copy information
*/ 
  for (i=0; i < NClusters; i++)
    {
      if ((Clusters[i].Son1 == -1) && (Clusters[i].Son2 == -1))
	{
/*
  Memory for the new blocks
*/
	  if (NBlocks == 0)
	    {
	      (*aFinBlocks)=(GCBlock *)malloc(1*sizeof(GCBlock));
	      assert((*aFinBlocks) != NULL);
	    }
	  else
	    {
	      (*aFinBlocks)=(GCBlock *)realloc((*aFinBlocks),(NBlocks+1)*sizeof(GCBlock));
	      assert((*aFinBlocks) != NULL);
	    }


	  (*aFinBlocks)[NBlocks].NBlock=Clusters[i].Number;
	  (*aFinBlocks)[NBlocks].MBlock=Blocks[NCheck].MBlock;
	  (*aFinBlocks)[NBlocks].Cls=Blocks[NCheck].Cls;
	  (*aFinBlocks)[NBlocks].Rows=(long *)malloc(Clusters[i].Number*sizeof(long));
	  assert((*aFinBlocks)[NBlocks].Rows != NULL);

	  for (j=0; j < Clusters[i].Number; j++)
	    {
	      if (Blocks[NCheck].Rows[Permu[Clusters[i].PermuPos+j]] > NGlobal)
		      printf("Alarm 2 in CLP : %ld\n",Blocks[NCheck].Rows[Permu[Clusters[i].PermuPos+j]]);
	      (*aFinBlocks)[NBlocks].Rows[j]=Blocks[NCheck].Rows[Permu[Clusters[i].PermuPos+j]];
	    }
	  (*aFinBlocks)[NBlocks].Cols=Blocks[NCheck].Cols;
      
	  if (((*aFinBlocks)[NBlocks].NBlock <= ParamACA.NCluMin) || ((*aFinBlocks)[NBlocks].MBlock <= ParamACA.NCluMin))
	    (*aFinBlocks)[NBlocks].Type=Dense;
	  else 
	    (*aFinBlocks)[NBlocks].Type=Admissible;

	  NBlocks+=1;
	}
    }
/*
  New number of final blocks
*/
  *aNBlocks=NBlocks;

  free(Clusters); free(Permu); free(X);
}
/*
 Cluster devision, unit sphere
*/
void GDevideCluster(double *X,long *aNClusters,long IClu,Cluster **aClusters, 
                    long *Permu,double Crit,ACAParameters ParamACA)
{
/*
  Local variables
*/
    long i,j,NClu,ClBegin,ClEnd,Num1,Num2,ICheck,IEnd,NClusters,IClu1,IClu2;
    double Mom[10],XMin[3],XMax[3],Centre[3],Diag[3],Diff[3],XT[3];
    double DiagLength,LCrit,alpha,Radius,Phi;
    Cluster *Clusters;
/*
  LAPACK variables
*/
    int N,LDA,LWORK,INFO;
    double A[9],W[3],WORK[8];
/*
  Number of elements
*/
    Clusters=*aClusters;
    NClusters=*aNClusters;
/*
  Number of elements
*/
    NClu=Clusters[IClu].Number;

    if (NClu <= ParamACA.NCluMin)
/*
  Small cluster, no devision
*/
    {
	Clusters[IClu].Son1=-1;
	Clusters[IClu].Son2=-1;
	return;
    }

    ClBegin=Clusters[IClu].PermuPos;
    ClEnd=Clusters[IClu].PermuPos+NClu-1;
/*
  Moments
*/
    dset(10,0.0,Mom,1);

    for(i=ClBegin; i <= ClEnd; i++)
      {
	j=Permu[i];

        Mom[0]+=1.0;
        Mom[1]+=X[3*j];
        Mom[2]+=X[3*j+1];
        Mom[3]+=X[3*j+2];
        Mom[4]+=X[3*j]*X[3*j];
        Mom[5]+=X[3*j]*X[3*j+1];
        Mom[6]+=X[3*j]*X[3*j+2];
        Mom[7]+=X[3*j+1]*X[3*j+1];
        Mom[8]+=X[3*j+1]*X[3*j+2];
        Mom[9]+=X[3*j+2]*X[3*j+2];
      }
/*
  Centre
*/
    cblas_dcopy(3,Mom+1,1,Centre,1);
    alpha=1.0/Mom[0];
    cblas_dscal(3,alpha,Centre,1);
    alpha=1.0/cblas_dnrm2(3,Centre,1);
    cblas_dscal(3,alpha,Centre,1);
/*
  Maximal angle
*/
    Phi=0.0;
   
    for(i=ClBegin; i <= ClEnd; i++)
      Phi=max(Phi,acos(cblas_ddot(3,X+3*Permu[i],1,Centre,1)));
/*
  Angle condition
*/
    if (Phi <= Crit)
/*
  Admissible cluster, no devision
*/
    {
 	Clusters[IClu].Son1=-1;
	Clusters[IClu].Son2=-1;
	return;
    }
/*
  Devision forced, covariance matrix
*/
    A[0]=Mom[4]-Centre[0]*Mom[1];
    A[1]=A[3]=Mom[5]-Centre[0]*Mom[2];
    A[2]=A[6]=Mom[6]-Centre[0]*Mom[3];
    A[4]=Mom[7]-Centre[1]*Mom[2];
    A[5]=A[7]=Mom[8]-Centre[1]*Mom[3];
    A[8]=Mom[9]-Centre[2]*Mom[3];
/*
  LAPACK preparations
*/
    N=3;
    LDA=3;
    LWORK=8;
/*
  LAPACK for eigevalues and eigenvectors
*/
    dsyev_("V","U",&N,A,&LDA,W,WORK,&LWORK,&INFO);
    if (INFO != 0)
    {
	printf("LAPACK Error in dsyev: INFO = %d\n",INFO);
        exit(0);
    }
/*
  Devision
*/
	LCrit=cblas_ddot(3,Centre,1,A+6,1);

        Num1=0;
        Num2=0;
 
        ICheck=ClBegin;
        IEnd=ClEnd;

        for (i=0; i < NClu; i++)
	{
            j=Permu[ICheck];

	    if (cblas_ddot(3,X+3*j,1,A+6,1) >= LCrit)
	    {
		Num1+=1;
                ICheck+=1;
	    } else {
                Num2+=1;
                Permu[ICheck]=Permu[IEnd];
		Permu[IEnd]=j;
		IEnd-=1;
	    }
	}
/*
  Save the first son
*/
        if ((Num1 == 0) || (Num2 == 0))
	{
/*
  Non-separable cluster
*/
	    Clusters[IClu].Son1=-1;
	    Clusters[IClu].Son2=-1;
	    printf("The grid contains non-separable objects !\n");
	    return;
	} else {
/*
  New clusters
*/
            NClusters+=2;
            
	    Clusters=realloc(Clusters,NClusters*sizeof(Cluster));
	    assert(Clusters != NULL);
/*
  Information and further devision
*/
	    IClu1=NClusters-2;
	    IClu2=NClusters-1;

            Clusters[IClu].Son1=IClu1;

            Clusters[IClu1].Level=Clusters[IClu].Level+1;            
            Clusters[IClu1].Father=IClu;            
            Clusters[IClu1].Number=Num1;            
            Clusters[IClu1].PermuPos=ClBegin;            

            Clusters[IClu].Son2=IClu2;
         
            Clusters[IClu2].Level=Clusters[IClu].Level+1;            
            Clusters[IClu2].Father=IClu;            
            Clusters[IClu2].Number=Num2;            
            Clusters[IClu2].PermuPos=ClBegin+Num1;            
	} 

	*aClusters=Clusters;
	*aNClusters=NClusters;
}
void CheckCand(double *Points,Cluster *Clusters,GCBlock *Blocks,
               long NCheck,ACAParameters ParamACA)
{
/*
  Local variables
*/
  long IEnd,RowCk,ih,ii;
  double R,Crit,Dist[3],C[3],DistNorm;
/*
  Preparations
*/
  IEnd=Blocks[NCheck].NBlock-1;
  RowCk=0;

  R=Clusters[Blocks[NCheck].Cls].Radius;
  Crit=ParamACA.ACAAlpha*R+ParamACA.ACABeta*ParamACA.ACAKappa*R*R;

  cblas_dcopy(3,Clusters[Blocks[NCheck].Cls].Centre,1,C,1);  

  do
    {
      cblas_dcopy(3,C,1,Dist,1);
      cblas_daxpy(3,-1.0,Points+3*Blocks[NCheck].Rows[RowCk],1,Dist,1);

      DistNorm=cblas_dnrm2(3,Dist,1);

      if(DistNorm >= Crit)
	{
	  RowCk+=1;
	}
      else
	{
	  ih=Blocks[NCheck].Rows[IEnd];
          Blocks[NCheck].Rows[IEnd]=Blocks[NCheck].Rows[RowCk];
          Blocks[NCheck].Rows[RowCk]=ih;
          IEnd-=1;
	}
    }
  while (RowCk <= IEnd);

  Blocks[NCheck].NBlock=RowCk;
}
