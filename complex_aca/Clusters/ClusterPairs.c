/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                27.12.2008                                   *
\*****************************************************************************/
/*
  Last change 07.01.2011
*/

/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

long ClusterPairs(Cluster *ClustersRow,Cluster *ClustersColumn,Pair **aPairs,long *aNMax,long *aMMax,ACAParameters ParamACA)
{
/*
  Local variables
*/
  long NCand,NCheck,ICheck1,ICheck2,NBlock,MBlock,NMax,MMax;
  Pair *Pairs;
/*
  Initial settings
*/
  NCand=1;
  NCheck=0;


  ICheck1=0;
  ICheck2=0;

  Pairs=(Pair *)malloc(NCand*sizeof(Pair));
  assert(Pairs != NULL);

  Pairs[0].Clu1=0;
  Pairs[0].Clu2=0;

  NMax=0;
  MMax=0;
/*
  Check all possible pairs
*/
  do
    {
      NBlock=ClustersRow[ICheck1].Number;
      MBlock=ClustersColumn[ICheck2].Number;
/*
  Dense block
*/
      if ((ClustersRow[ICheck1].Son1 == -1) || (ClustersColumn[ICheck2].Son1 == -1))
	{
	  Pairs[NCheck].Type=Dense;

	  NMax=lmax(NBlock,NMax);
	  MMax=lmax(MBlock,MMax);
	}
/*
  Admissible block
*/
      else if(AdmissiblePair(ICheck1,ICheck2,ClustersRow,ClustersColumn,ParamACA) == 1)
	{
	  Pairs[NCheck].Type=Admissible;

	  NMax=lmax(NBlock,NMax);
	  MMax=lmax(MBlock,MMax);
	}
/*
  Hierarchical block
*/
      else 
	{
	  Pairs[NCheck].Type=Hierarchical;

	  Pairs[NCheck].A11=NCand;
	  Pairs[NCheck].A12=NCand+1;
	  Pairs[NCheck].A21=NCand+2;
	  Pairs[NCheck].A22=NCand+3;

	  Pairs=(Pair *)realloc(Pairs,(NCand+4)*sizeof(Pair));
	  assert(Pairs != NULL);

	  Pairs[NCand].Clu1=ClustersRow[ICheck1].Son1;
	  Pairs[NCand].Clu2=ClustersColumn[ICheck2].Son1;

	  Pairs[NCand+1].Clu1=ClustersRow[ICheck1].Son1;
	  Pairs[NCand+1].Clu2=ClustersColumn[ICheck2].Son2;

	  Pairs[NCand+2].Clu1=ClustersRow[ICheck1].Son2;
	  Pairs[NCand+2].Clu2=ClustersColumn[ICheck2].Son1;

	  Pairs[NCand+3].Clu1=ClustersRow[ICheck1].Son2;
	  Pairs[NCand+3].Clu2=ClustersColumn[ICheck2].Son2;

	  NCand+=4;

/* 	  Pairs=(Pair *)realloc(Pairs,(NCand+2)*sizeof(Pair)); */
/* 	  assert(Pairs != NULL); */

/* 	  if(NBlock >= MBlock) */
/* 	    { */
/* 	      Pairs[NCand].Clu1=ClustersRow[ICheck1].Son1; */
/* 	      Pairs[NCand].Clu2=ICheck2; */
    
/* 	      Pairs[NCand+1].Clu1=ClustersRow[ICheck1].Son2; */
/* 	      Pairs[NCand+1].Clu2=ICheck2; */
/* 	    }  */
/* 	  else  */
/* 	    { */
/* 	      Pairs[NCand].Clu1=ICheck1; */
/* 	      Pairs[NCand].Clu2=ClustersColumn[ICheck2].Son1; */
    
/* 	      Pairs[NCand+1].Clu1=ICheck1; */
/* 	      Pairs[NCand+1].Clu2=ClustersColumn[ICheck2].Son2; */
/* 	    } */
/* 	    NCand+=2;          */
	}

      NCheck+=1;

      if (NCheck < NCand)
	{
	  ICheck1=Pairs[NCheck].Clu1;
	  ICheck2=Pairs[NCheck].Clu2;
	}
    }
  while (NCheck < NCand);

  *aPairs=Pairs;
  *aNMax=NMax;
  *aMMax=MMax;

  return NCand;    
}
