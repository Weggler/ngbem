/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                  Adaptive Cross Approximation for EMSS                      *
*                                                                             *
*                                16.02.2011                                   *
\*****************************************************************************/
/*
  Last change 11.03.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void BHMatrixToBlock(BlockHMatrix *BHMat,CBlock **A)
{
  /* Local variables */
  long NB,MB,N,M;
  long N2,ND,i,j,ip,ib,iab;

  CBlock **P;
  CBlock *AB;

  /* Block dimension of the matrix */
  NB=(*BHMat).NBlockRow; 
  MB=(*BHMat).NBlockColumn;
  assert(NB == MB);

  /* Dimension of the matrix */
  N=(*BHMat).NRow; 
  M=(*BHMat).NColumn;
  assert(N == M);

  if (NB == 1)
  {
    *A=((*BHMat).HBlocks[0]).CBlocks;
  }
  else
  {
    /* Memory allocation for the A Block */
    *A=(CBlock *)malloc(1*sizeof(CBlock));
    assert(*A != NULL);
    
    /* Next power of 2 */
    N2=(long)pow(2,(long)(log(NB-0.5)/log(2.0)+1));
//    printf("BHMatrixToBlock: N2 = %d, NB = %d\n", N2, NB);
    
    /* Memory for the pointers */
    P=(CBlock **)malloc(N2*N2*sizeof(CBlock *));
    assert(P != NULL);
//    printf("BHMatrixToBlock: N2*N2 = %d\n", N2*N2);
    
    /* Additional blocks */
    ND=N2-NB;
    if (ND > 0)
    {
      AB=(CBlock *)malloc(ND*(MB+NB+ND)*sizeof(CBlock));
      assert(AB != NULL);
    }      
    
    /* Copy the pointers */
    ip =0;
    ib =0;
    iab=0;
    
    for (j=0; j < MB; j++)
    {
      for (i=0; i < NB; i++)
      {
        P[ip]=((*BHMat).HBlocks[ib]).CBlocks;
//        printf("BHMatrixToBlock: %d, %d, %d\n", ip,ib,((*BHMat).HBlocks[ib]).NBlocks);
//        fflush(stdout);
//        printf("BHMatrixToBlock: %p\n", ((*BHMat).HBlocks[ib]).CBlocks);
//        printf("BHMatrixToBlock: %ld\n", (*P[ip]).IBlock);
//        fflush(stdout);
        ip++;
        ib++;
      }   
      /* Additional zero rows */
      if (ND > 0)
      {
        for (i=0; i < ND; i++)
        {
          AB[iab].IBlock=N+i;
          AB[iab].JBlock=(*BHMat).MDim[j];
          AB[iab].NBlock=1;
          if (j < MB-1)
          {
            AB[iab].MBlock=(*BHMat).MDim[j+1]-(*BHMat).MDim[j];
//            printf("BHMatrixToBlock 1: AB[%ld].MBlock = %ld\n", iab, AB[iab].MBlock);
//            printf(".............from  %ld - %ld\n",(*BHMat).MDim[j+1],(*BHMat).MDim[j]);
          }
          else
          {
            AB[iab].MBlock=M-(*BHMat).MDim[j];
//            printf("BHMatrixToBlock 2: AB[%ld].MBlock = %ld\n", iab, AB[iab].MBlock);
//            printf(".............from  %ld - %ld\n",M,(*BHMat).MDim[j]);
          }
          
          AB[iab].Type=Null;
          
          AB[iab].A11=NULL;
          AB[iab].A21=NULL;
          AB[iab].A12=NULL;
          AB[iab].A22=NULL;
          
          AB[iab].Data=NULL;
          
          AB[iab].Rank=0;
          
          AB[iab].UData=NULL;
          AB[iab].VData=NULL;
          
          P[ip]=AB+iab;
          ip ++;
          iab++;
      	}   
      }
    }
    /* Additional zero columns */
    if (ND > 0)
    {
      for (j=0; j < ND; j++)
      {
        for (i=0; i < NB; i++)
        {
          AB[iab].IBlock=(*BHMat).NDim[i];
          AB[iab].JBlock=M+j;
          if (i < NB-1)
          {
            AB[iab].NBlock=(*BHMat).NDim[i+1]-(*BHMat).NDim[i];
//            printf("BHMatrixToBlock 3: AB[%ld].NBlock = %ld\n", iab, AB[iab].NBlock);
//            printf(".............from  %ld - %ld\n",(*BHMat).NDim[i+1],(*BHMat).NDim[i]);
          }
          else
          {
            AB[iab].NBlock=N-(*BHMat).NDim[i];
//            printf("BHMatrixToBlock 4: AB[%ld].NBlock = %ld\n", iab, AB[iab].NBlock);
//            printf(".............from  %ld - %ld\n",N,(*BHMat).NDim[i]);
          }
          AB[iab].MBlock=1;
        
          AB[iab].Type=Null;
        
          AB[iab].A11=NULL;
          AB[iab].A21=NULL;
          AB[iab].A12=NULL;
          AB[iab].A22=NULL;
        
          AB[iab].Data=NULL;
        
          AB[iab].Rank=0;
        
          AB[iab].UData=NULL;
          AB[iab].VData=NULL;
        
          P[ip]=AB+iab;
          ip++;
          iab++;
        }
        /* Additional 1x1 blocks */
        for (i=0; i < ND; i++)
        {
           AB[iab].IBlock=N+i;
           AB[iab].JBlock=M+j;
           AB[iab].NBlock=1;
           AB[iab].MBlock=1;
           
           AB[iab].A11=NULL;
           AB[iab].A21=NULL;
           AB[iab].A12=NULL;
           AB[iab].A22=NULL;
           
           if (i != j)
           {
             AB[iab].Type=Null;
             AB[iab].Data=NULL;
           }
           else
           {
             AB[iab].Type=Dense;
             AB[iab].Data=(double complex *)malloc(1*sizeof(double complex));
             assert(AB[iab].Data != NULL);
             AB[iab].Data[0]=1.0+I*0.0;
           }
           
           AB[iab].Rank=0;
           
           AB[iab].UData=NULL;
           AB[iab].VData=NULL;
           
           P[ip]=AB+iab;
           ip++;
           iab++;
        }
      }
    }
    PointerListToBlock(0,0,N+ND,M+ND,N2,P,*A);
    /* Free memory */
    free(P);
  }
}
