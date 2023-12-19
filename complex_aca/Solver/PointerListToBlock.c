/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                16.02.2011                                   *
\*****************************************************************************/
/*
  Last change 16.02.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void PointerListToBlock(long IB, long JB, 
                        long N,  long M,  long NB, 
                        CBlock **P, 
                        CBlock  *A)
{
/*
  Local variables
*/
  long NH;
  long N11,M11,N21,M21,N12,M12,N22,M22;
  long i,j,ip,ipa;

  CBlock **PA;
//  printf("Block data : %ld %ld %ld %ld %ld\n",IB,JB,N,M,NB); 
//  for (j=0; j < NB; j++) 
//  {
//    for (i=0; i < NB; i++) 
//    { 
//      printf(" Block : %ld %ld\n",i,j); 
//      printf("      Type              : %ld\n",(long)(*P[j*NB+i]).Type); 
//      printf("      Position          : %ld %ld\n",(*P[j*NB+i]).IBlock,(*P[j*NB+i]).JBlock); 
//      printf("      Dimension         : %ld %ld\n",(*P[j*NB+i]).NBlock,(*P[j*NB+i]).MBlock); 
//      printf("      Blocks            : %p %p %p %p\n",(*P[j*NB+i]).A11,(*P[j*NB+i]).A21,(*P[j*NB+i]).A12,(*P[j*NB+i]).A22); 
//      printf("      Data              : %p\n",(*P[j*NB+i]).Data); 
//      printf("      Rank              : %ld\n",(*P[j*NB+i]).Rank); 
//      printf("      UData             : %p\n",(*P[j*NB+i]).UData); 
//      printf("      VData             : %p\n",(*P[j*NB+i]).VData); 
//    } 
//  } 
//  fflush(stdout); 
//  getchar(); 

  (*A).IBlock = IB;
  (*A).JBlock = JB;
  (*A).NBlock = N;
  (*A).MBlock = M;
  (*A).Type   = Hierarchical;
  (*A).Data   = NULL;
  (*A).Rank   = 0;
  (*A).UData  = NULL;
  (*A).VData  = NULL;

/*
  Blocks
*/
  if(NB == 2)
  {
    (*A).A11=P[0];
    (*A).A21=P[1];
    (*A).A12=P[2];
    (*A).A22=P[3];
  }
  else
  {
    NH=NB/2;
/*
  Memory for the additional pointers
*/
    PA=(CBlock **)malloc(NH*NH*sizeof(CBlock *));
    assert(PA != NULL);
/*
  Block 1,1
*/
    (*A).A11=(CBlock *)malloc(1*sizeof(CBlock));
    assert((*A).A11 != NULL);
  
    N11=0;
    M11=0;   
	  
    for (i=0; i < NH; i++)
    {
      N11+=(*P[i]).NBlock;
    }
    
    for (j=0; j < NH; j++)
    {
      M11+=(*P[j*NB]).MBlock;
    }
     
    ipa=0;
    ip =0; 
     
    for (j=0; j < NH; j++)
    {
      for (i=0; i < NH; i++)
      {
        PA[ipa]=P[ip];
        ipa++;
        ip++;
      }
      ip+=NH;
    }

    PointerListToBlock(IB,JB,N11,M11,NH,PA,(*A).A11);
/*
  Block 2,1
*/
    (*A).A21=(CBlock *)malloc(1*sizeof(CBlock));
    assert((*A).A21 != NULL);
    
    N21=N-N11;
    M21=M11;
    
    ipa=0;
    ip=NH;     
    
    for (j=0; j < NH; j++)
    {
      for (i=NH; i < NB; i++)
      {
        PA[ipa]=P[ip];
        ipa++;
        ip++;
      }
      ip+=NH;
    }

    PointerListToBlock(IB+N11,JB,N21,M21,NH,PA,(*A).A21);
/*
  Block 1,2
*/
    (*A).A12=(CBlock *)malloc(1*sizeof(CBlock));
    assert((*A).A12 != NULL);
    
    N12=N11;
    M12=M-M11;
    
    ipa=0;
    ip =NB*NH;     
    
    for (j=NH; j < NB; j++)
    {
      for (i=0; i < NH; i++)
      {
        PA[ipa]=P[ip];
        ipa++;
        ip++;
      }
      ip+=NH;
    }

    PointerListToBlock(IB,JB+M11,N12,M12,NH,PA,(*A).A12);
/*
  Block 2,2
*/
    (*A).A22=(CBlock *)malloc(1*sizeof(CBlock));
    assert((*A).A22 != NULL);
  
    N22=N-N11;
    M22=M-M11;
  
    ipa=0;
    ip =NB*NH+NH;     
  
    for (j=NH; j < NB; j++)
    {
      for (i=NH; i < NB; i++)
      {
        PA[ipa]=P[ip];
        ipa++;
        ip++;
      }
      ip+=NH;
    }

    PointerListToBlock(IB+N11,JB+M11,N22,M22,NH,PA,(*A).A22);
/*
  Free memory
*/
    free(PA);
  }
}

