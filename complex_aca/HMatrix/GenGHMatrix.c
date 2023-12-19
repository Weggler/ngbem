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

void GenGHMatrix(long N,long M,
                 long NBlocks,GCBlock *Blocks,GHMatrix *aGHMat,
                 double complex (*ElementName)(long,long),
                 ACAParameters *aParamACA)
{
/*
  Local variables
*/
  long IBlock,NBlock,MBlock,Result,Rank,NBlDense,NBlACA,NBlFailed;
  double complex *EBlock,*UBlock,*VBlock;
  GHMatrix GHMat;
  BlockType BlType;
  ACAParameters ParamACA;

  long const True=1,False=0;
/*
  Initialisation
*/
  NBlDense=0;
  NBlACA=0;
  NBlFailed=0;

  ParamACA=*aParamACA;
  GHMat=*aGHMat;

  GHMat.NRow=N; GHMat.NColumn=M;
  GHMat.NBlocks=NBlocks;
  GHMat.GCBlocks=Blocks;
/*
  Generate blocks
*/
  ParamACA.Memory=0.0;
  ParamACA.AppMemory=0.0;

  for (IBlock=0; IBlock < NBlocks; IBlock++)
    {
      BlType=Blocks[IBlock].Type;
/*
      printf("%ld %ld %ld\n",IBlock,Blocks[IBlock].NBlock,Blocks[IBlock].MBlock);
*/
      NBlock=Blocks[IBlock].NBlock;
      MBlock=Blocks[IBlock].MBlock;
     
      if (BlType == Dense)
	{
          NBlDense+=1;
/*
  Memory without approximation
*/
	  ParamACA.Memory+=(double)NBlock*(double)MBlock*16.0/(1024.0*1024.0);
/*
  Memory allocation for dense block
*/
	  EBlock=(double complex *)malloc(NBlock*MBlock*sizeof(double complex));
	  assert(EBlock != NULL);
/*
  Generate dense block
*/
	  GenGBlock(Blocks[IBlock],0,0,NBlock,MBlock,EBlock,ElementName);

	  Blocks[IBlock].Data=EBlock;
   
	  ParamACA.AppMemory+=(double)NBlock*(double)MBlock*16.0/(1024.0*1024.0);
	}
      else if (BlType == Admissible)
	{
/*
  Memory without approximation
*/
	  ParamACA.Memory+=(double)NBlock*(double)MBlock*16.0/(1024.0*1024.0);
/*
  Memory allocation for low rank block
*/
	  UBlock=(double complex *)malloc(NBlock*ParamACA.MaxRank*sizeof(double complex));
	  assert(UBlock != NULL);

	  VBlock=(double complex *)malloc(MBlock*ParamACA.MaxRank*sizeof(double complex));
	  assert(VBlock != NULL);
/*
  Generate low rank matrix
*/
	  PartialGACA(Blocks[IBlock],
		      &Result,
		      &Rank,UBlock,VBlock,
		      ElementName,ParamACA);         
/*
  Memory reallocation
*/
	  if (Result == True)
	    {
             NBlACA+=1;
/*
  Post compression
*/

/*               printf("Old rank = %ld\n",Rank); */
/*               fflush(stdout);  */
              PostCompression(NBlock,MBlock,&Rank,UBlock,VBlock,ParamACA.ACAEps);
/*               printf("New rank = %ld\n",Rank); */
/*               getchar();  */
/*
  Memory reaalocation
*/
	      UBlock=(double complex *)realloc(UBlock,NBlock*Rank*sizeof(double complex));
	      assert(UBlock != NULL);

	      VBlock=(double complex *)realloc(VBlock,MBlock*Rank*sizeof(double complex));
	      assert(VBlock != NULL);

	      Blocks[IBlock].Rank=Rank;
	      Blocks[IBlock].UData=UBlock;
	      Blocks[IBlock].VData=VBlock;
   
	      ParamACA.AppMemory+=(double)Rank*((double)NBlock+(double)MBlock)*16.0/(1024.0*1024.0);
	    }
	  else
	    {
      printf("%ld %ld %ld\n",IBlock,Blocks[IBlock].NBlock,Blocks[IBlock].MBlock);
/*
       getchar();
*/
              NBlFailed+=1;

	      free(UBlock); free(VBlock);

	      Blocks[IBlock].Type=Dense;

	      EBlock=(double complex *)malloc(NBlock*MBlock*sizeof(double complex));
	      assert(EBlock != NULL);
/*
  Generate dense block
*/
	      GenGBlock(Blocks[IBlock],0,0,NBlock,MBlock,EBlock,ElementName);

	      Blocks[IBlock].Data=EBlock;
  
	      ParamACA.AppMemory+=(double)NBlock*(double)MBlock*16.0/(1024.0*1024.0);
	    }
	}
/*
  Progress
*/
      if((ParamACA.Info == 1) && (ParamACA.Memory > 0.0))
	{
	  printf("  %6.2f %% done  %8.2f MB instead of %8.2f MB (%6.2f %%)\n",
		 ParamACA.Memory/
		 (((double)N)*((double)M)*16.0/(1024.0*1024.0))*100.0,
		 ParamACA.AppMemory,ParamACA.Memory,
		 ParamACA.AppMemory/ParamACA.Memory*100.0);
	  fflush(stdout);
	}
    }

*aGHMat=GHMat;
*aParamACA=ParamACA;

printf("Blocks: Dense = %ld, ACA = %ld, Failed = %ld, All = %ld\n",NBlDense,NBlACA,NBlFailed,NBlDense+NBlACA+NBlFailed);
}
