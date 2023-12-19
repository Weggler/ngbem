//
//   S. Rjasanow: Adaptive Cross Approximation
//
//----------------------------------------------------------------------
//
//   routine name        - SymmClusterPairs
//
//----------------------------------------------------------------------
//
//   latest revision    - Aug 12
//
//   purpose            - compute symmetric cluster pairs for ACA
//
//   in:
//             Clusters - i.g. list of cluster trees --- 
//                        here we need only the one related to
//                        elements (==mdle nodes)
//
//   out
//               aPairs - list of cluster Pairs derived from cluster 
//                        trees
//                aNMax - list of max    row-dim wrt all (?) pairs
//                aMMax - list of max column-dim wrt all (?) pairs
//             ParamACA - all ACA parameters 
//
//----------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

long SymmClusterPairs(Cluster *Clusters,
                      Pair **aPairs,long *aNMax,long *aMMax,ACAParameters ParamACA)
{
  /* Local variables */
  long NCand,NCheck,ICheck1,ICheck2,NBlock,MBlock,NMax,MMax;
  Pair *Pairs;

  /* Initial settings */
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

  /* Check all possible pairs */
  do
  {
    NBlock=Clusters[ICheck1].Number;
    MBlock=Clusters[ICheck2].Number;
    
    /* Small block */
    if ((Clusters[ICheck1].Son1 == -1) || (Clusters[ICheck2].Son1 == -1))
    {
      Pairs[NCheck].Type=Dense;
    
      NMax=lmax(NBlock,NMax);
      MMax=lmax(MBlock,MMax);
    }

    /* Admissible block */
    else if (AdmissiblePair(ICheck1,ICheck2,Clusters,Clusters,ParamACA) == 1)
    {
      Pairs[NCheck].Type=Admissible;
    
      NMax=lmax(NBlock,NMax);
      MMax=lmax(MBlock,MMax);
    } 

    /* Hierarchical block */
    else 
    {
      Pairs[NCheck].Type=Hierarchical;
    
      Pairs[NCheck].A11=NCand;
      Pairs[NCheck].A12=NCand+1;
      Pairs[NCheck].A21=NCand+2;
      Pairs[NCheck].A22=NCand+3;
    
      Pairs=(Pair *)realloc(Pairs,(NCand+4)*sizeof(Pair));
      assert(Pairs != NULL);
    
      Pairs[NCand].Clu1=Clusters[ICheck1].Son1;
      Pairs[NCand].Clu2=Clusters[ICheck2].Son1;
    
      Pairs[NCand+1].Clu1=Clusters[ICheck1].Son1;
      Pairs[NCand+1].Clu2=Clusters[ICheck2].Son2;
    
      Pairs[NCand+2].Clu1=Clusters[ICheck1].Son2;
      Pairs[NCand+2].Clu2=Clusters[ICheck2].Son1;
    
      Pairs[NCand+3].Clu1=Clusters[ICheck1].Son2;
      Pairs[NCand+3].Clu2=Clusters[ICheck2].Son2;
    
      NCand+=4;
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
