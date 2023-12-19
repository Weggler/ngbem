//------------------------------------------------------------------
//
//  routine name  - GenBlock
//
//------------------------------------------------------------------
//
//  last revision - Aug 12
//
//  purpose       - generate block hierarchic matrix
//
//  in:
//  out:
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

void GenBlock(Pair BlockPair,
              long IBPos,long JBPos,
              long *PermuRow,long *PermuColumn,
              Cluster *ClustersRow, Cluster *ClustersColumn,
              long IBlock,long JBlock,
              long NBlock,long MBlock,
              double complex *Block,
              double complex (*ElementName)(long,long))
{
/* Local variables */
  long IB,JB,IEl,JEl,IBegin,JBegin;

/* Initialisation */
  IBegin=ClustersRow[BlockPair.Clu1].PermuPos; //!= IBlock from Input
  JBegin=ClustersColumn[BlockPair.Clu2].PermuPos; //!= JBlock

/* Generation columnwise */
  for (JB=JBlock; JB < JBlock+MBlock; JB++)
  {
//    JEl=PermuColumn[JBegin+JB]+JBPos;
    JEl=PermuColumn[JBegin+JB];
 
    for (IB=IBlock; IB < IBlock+NBlock; IB++)
    {
//      IEl=PermuRow[IBegin+IB]+IBPos;
      IEl=PermuRow[IBegin+IB];
 
      Block[NBlock*(JB-JBlock)+IB-IBlock]=ElementName(IEl,JEl);
    }      
  }
}
