//
//  routine name  -  SaveHMatrix
//
//------------------------------------------------------------------
//
//  last revision -  Aug 12
//
//  purpose       -  writes a complex symmetric HMatrix in a file
//
//  in                                                              
//           HMat - the HMatrix to be saved
//HMatrixFileName - file name
//
//  out           - HMatrixFileName filled with data 
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/GlobalVariables.h"

void SaveCSymHMatrix(HMatrix HMat,char *HMatrixFileName)
{
  /*
     Local variables
   */
  int Ll,Ldc,Lpt,Lbt,FileD,WriteD;
  FILE *HMatrixFile;
  long IPair;
  /*
     Initialisation
   */
  Ll=sizeof(long);
  Ldc=sizeof(double complex);
  Lpt=sizeof(PrecondType);
  Lbt=sizeof(BlockType);
  /*
     Open file
  */
  FileD=open(HMatrixFileName,O_WRONLY | O_CREAT, S_IRWXU);
//  FileD=open(HMatrixFileName,O_WRONLY | O_CREAT);
  assert(FileD != -1);
  /*
     Write dimensions
  */
  WriteD=write(FileD,&HMat.NRow,Ll);
  assert(WriteD != -1);
  
  WriteD=write(FileD,&HMat.NColumn,Ll);
  assert(WriteD != -1);
  /*
     Write permutations
  */
  WriteD=write(FileD,HMat.PermuRow,HMat.NRow*Ll);
//  fprintf(stderr,"Save NRow= %ld\n",HMat.NRow);
  assert(WriteD != -1);
  
  WriteD=write(FileD,HMat.PermuColumn,HMat.NColumn*Ll);
  assert(WriteD != -1);
  /*
     Write number of blocks
  */
  WriteD=write(FileD,&HMat.NBlocks,Ll);
  assert(WriteD != -1);
  /*
     Write type of preconditioning
  */
  WriteD=write(FileD,&HMat.Precond,Lpt);
  assert(WriteD != -1);
  /*
     Write blocks
  */
  for (IPair=0; IPair < HMat.NBlocks; IPair++)
  {
    /*
       Write position
    */
    WriteD=write(FileD,&(HMat.CBlocks[IPair].IBlock),Ll);
    assert(WriteD != -1);
    
    WriteD=write(FileD,&(HMat.CBlocks[IPair].JBlock),Ll);
    assert(WriteD != -1);
    /*
       Write dimension
    */
    WriteD=write(FileD,&(HMat.CBlocks[IPair].NBlock),Ll);
    assert(WriteD != -1);
    
    WriteD=write(FileD,&(HMat.CBlocks[IPair].MBlock),Ll);
    assert(WriteD != -1);
    /*
       Write type
    */
    WriteD=write(FileD,&(HMat.CBlocks[IPair].Type),Lbt);
    assert(WriteD != -1);
    
    long bn;
    switch(HMat.CBlocks[IPair].Type)
    {
      case Hierarchical:
        bn=HMat.CBlocks[IPair].A11 - HMat.CBlocks; 
        WriteD=write(FileD,&bn,Ll);
        assert(&(HMat.CBlocks[bn])==HMat.CBlocks[IPair].A11);
        bn=HMat.CBlocks[IPair].A12 - HMat.CBlocks; 
        WriteD=write(FileD,&bn,Ll);
        assert(&(HMat.CBlocks[bn])==HMat.CBlocks[IPair].A12);
        bn=HMat.CBlocks[IPair].A21 - HMat.CBlocks; 
        WriteD=write(FileD,&bn,Ll);
        assert(&(HMat.CBlocks[bn])==HMat.CBlocks[IPair].A21);
        bn=HMat.CBlocks[IPair].A22 - HMat.CBlocks; 
        WriteD=write(FileD,&bn,Ll);
        assert(&(HMat.CBlocks[bn])==HMat.CBlocks[IPair].A22);
        break;

      case Dense:
        if(HMat.CBlocks[IPair].IBlock <= HMat.CBlocks[IPair].JBlock)
        {
          WriteD=write(FileD,HMat.CBlocks[IPair].Data,
                 HMat.CBlocks[IPair].NBlock*HMat.CBlocks[IPair].MBlock*Ldc);
          assert(WriteD != -1);
        }
        break;

      case Admissible:
        WriteD=write(FileD,&(HMat.CBlocks[IPair].Rank),Ll);
        assert(WriteD != -1);
        if(HMat.CBlocks[IPair].IBlock <= HMat.CBlocks[IPair].JBlock)
        {
          WriteD=write(FileD,HMat.CBlocks[IPair].UData,
                 HMat.CBlocks[IPair].NBlock*HMat.CBlocks[IPair].Rank*Ldc);
          assert(WriteD != -1);

          WriteD=write(FileD,HMat.CBlocks[IPair].VData,
                 HMat.CBlocks[IPair].MBlock*HMat.CBlocks[IPair].Rank*Ldc);
          assert(WriteD != -1);
        }
        break;
    }
  }
  close(FileD);
}

void LoadCSymHMatrix(HMatrix* HMat,char *HMatrixFileName)
{
  /*
     Local variables
  */
  int Ll,Ldc,Lpt,Lbt,FileD,ReadD;
  FILE *HMatrixFile;
  long IPair;
  /*
     Initialisation
  */
  Ll=sizeof(long);
  Ldc=sizeof(double complex);
  Lpt=sizeof(PrecondType);
  Lbt=sizeof(BlockType);
  /*
     Open file
  */
  FileD=open(HMatrixFileName,O_RDONLY);
  assert(FileD != -1);
  /*
     Read dimensions
  */
  ReadD=read(FileD,&HMat->NRow,Ll);
//  fprintf(stderr,"Load NRow= %ld\n",HMat->NRow);
  assert(ReadD != -1);
  
  ReadD=read(FileD,&HMat->NColumn,Ll);
  assert(ReadD != -1);
  /*
     Read permutations
  */
  HMat->PermuRow=malloc(HMat->NRow*Ll);
  ReadD=read(FileD,HMat->PermuRow,HMat->NRow*Ll);
  assert(ReadD != -1);
  
  HMat->PermuColumn=malloc(HMat->NColumn*Ll);
  ReadD=read(FileD,HMat->PermuColumn,HMat->NColumn*Ll);
  assert(ReadD != -1);
  /*
     Read number of blocks
  */
  ReadD=read(FileD,&HMat->NBlocks,Ll);
  assert(ReadD != -1);
  /*
     Read type of preconditioning
  */
  ReadD=read(FileD,&HMat->Precond,Lpt);
  assert(ReadD != -1);
  HMat->CBlocks=malloc( HMat->NBlocks *sizeof(CBlock));
  /*
     Read blocks
  */
  for (IPair=0; IPair < HMat->NBlocks; IPair++)
  {
    /*
       Read position
    */
    ReadD=read(FileD,&(HMat->CBlocks[IPair].IBlock),Ll);
    assert(ReadD != -1);
    
    ReadD=read(FileD,&(HMat->CBlocks[IPair].JBlock),Ll);
    assert(ReadD != -1);
    /*
       Read dimension
    */
    ReadD=read(FileD,&(HMat->CBlocks[IPair].NBlock),Ll);
    assert(ReadD != -1);
    
    ReadD=read(FileD,&(HMat->CBlocks[IPair].MBlock),Ll);
    assert(ReadD != -1);
    /*
       Read type
    */
    ReadD=read(FileD,&(HMat->CBlocks[IPair].Type),Lbt);
    assert(ReadD != -1);
    /*
       Complex Symmetric HMatrix
    */
    HMat->Symmetry = ComplexSymmetric;
    /*
       Conventional ACA
    */
    HMat->NLocColumn = NULL; 
    HMat->NLocRow = NULL; 

    long bn;
    switch(HMat->CBlocks[IPair].Type)
    {
      case Hierarchical:
        ReadD=read(FileD,&bn,Ll); 
        HMat->CBlocks[IPair].A11 = &(HMat->CBlocks[bn]);
        ReadD=read(FileD,&bn,Ll); 
        HMat->CBlocks[IPair].A12 = &(HMat->CBlocks[bn]);
        ReadD=read(FileD,&bn,Ll); 
        HMat->CBlocks[IPair].A21 = &(HMat->CBlocks[bn]);
        ReadD=read(FileD,&bn,Ll); 
        HMat->CBlocks[IPair].A22 = &(HMat->CBlocks[bn]);
        break;
      
      case Dense:
        if(HMat->CBlocks[IPair].IBlock <= HMat->CBlocks[IPair].JBlock)
        {
          HMat->CBlocks[IPair].Data=
          malloc(HMat->CBlocks[IPair].NBlock*HMat->CBlocks[IPair].MBlock*Ldc);
          ReadD=read(FileD,HMat->CBlocks[IPair].Data,
                HMat->CBlocks[IPair].NBlock*HMat->CBlocks[IPair].MBlock*Ldc);
          assert(ReadD != -1);
        }
        break;
      
      case Admissible:
        ReadD=read(FileD,&(HMat->CBlocks[IPair].Rank),Ll);
        assert(ReadD != -1);
        if(HMat->CBlocks[IPair].IBlock <= HMat->CBlocks[IPair].JBlock)
        {
          HMat->CBlocks[IPair].UData=
          malloc(HMat->CBlocks[IPair].NBlock*HMat->CBlocks[IPair].Rank*Ldc);
          ReadD=read(FileD,HMat->CBlocks[IPair].UData,
                HMat->CBlocks[IPair].NBlock*HMat->CBlocks[IPair].Rank*Ldc);
          assert(ReadD != -1);
      
          HMat->CBlocks[IPair].VData=
          malloc(HMat->CBlocks[IPair].MBlock*HMat->CBlocks[IPair].Rank*Ldc);
          ReadD=read(FileD,HMat->CBlocks[IPair].VData,
                HMat->CBlocks[IPair].MBlock*HMat->CBlocks[IPair].Rank*Ldc);
          assert(ReadD != -1);
        }
        break;
    }
  }
  close(FileD);
}
