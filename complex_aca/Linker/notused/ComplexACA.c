/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                05.05.2009                                   *
\*****************************************************************************/
/*
  Includes for the main function 
*/

# include "IncludesMain.h"

int main()
{
/*
  Local variables
*/
    char Directory[4096],Command[4096],ResultName[4096];
    long Error,MyId,NumThreads;
    FILE *ResultFile;
/*
  Project variables
*/
    Node *Nodes;
    Triangle *Triangles;
/* 
   Read input file 
*/
    Error=ReadInputFile();
    if (Error != 0)
    { 
	printf("Errors reading input file Input.dat !\n");
	exit(0);
    }
/* 
   Create output directory
*/
    sprintf(Directory,"%s%s/",Prefix,Suffix);

    sprintf(Command,"mkdir %s",Directory);
    system(Command);
/* 
   Create Info.dat file
*/
    sprintf(Command,"date > %sInfo.dat",Directory);
    system(Command);

    sprintf(Command,"echo >> %sInfo.dat",Directory);
    system(Command);
/* 
   Copy Input.dat in Info.dat file
*/
    sprintf(Command,"cat Input.dat >> %sInfo.dat",Directory);
    system(Command);

    sprintf(Command,"echo >> %sInfo.dat",Directory);
    system(Command);
/* 
   Allocation of the memory for the Nodes and Triangles
*/
    Nodes=(Node *)malloc(NumNodes*sizeof(Node));
    assert(Nodes != NULL);

    Triangles=(Triangle *)malloc(NumElements*sizeof(Triangle));
    assert(Triangles != NULL);
/*
  Simulation
*/
    Simulation(Nodes,Triangles);
/* 
   Complete Info.dat file
*/
    sprintf(Command,"date >> %sInfo.dat",Directory);
    system(Command);

    return 0;
}
/*****************************************************************************\
 *                                                                             *
 *                                  Simulation                                 *
 *                                                                             *
\*****************************************************************************/
void Simulation(Node *Nodes, Triangle *Triangles)
{
/*
  Local variables
*/
  ACAParameters ParamACA;
  long i,Error,NClusters,NPairs,IClu;
  Cluster *Clusters;
  Pair *Pairs;
  long *Permu,NMax,MMax,i1,i2,i3;
  double alpha,*X,*G;
  HMatrix HMat;
  BlockHMatrix BHMat;
  double complex *xv,*yv;
  FILE *ResultFile;
/* 
   Read geometry file
*/
  printf("Read geometry ...\n"); 
  fflush(stdout);

  ReadGeometry(Nodes,Triangles);

  printf("  Done with %ld Nodes and %ld Elements\n\n",NumNodes,NumElements);
  fflush(stdout);
/* 
  Construct cluster tree for nodes
*/
  G=(double *)malloc(NumElements*sizeof(double));
  assert(G != NULL);

  X=(double *)malloc(3*NumElements*sizeof(double));
  assert(X != NULL);
/*
   Find the neighbours of the nodes
*/
  printf("Find neighbours ...\n");
  fflush(stdout);

  FindNodeNeighbours(Nodes,Triangles);

  printf("  Done\n\n");
  fflush(stdout);
/*
  Compute geometrical information for the Triangles
*/
  printf("Compute geometrical data ...\n");
  fflush(stdout);

  Error=ElementInfo(Triangles);
  if (Error != 0)
    {
      printf("Error in mesh !\n");
      exit(0);
    }

  printf("  Done\n\n");
  fflush(stdout);
/*
   Prepare cluster data for elements
*/
  for (i=0; i < NumElements; i++)
    {
      i1=Triangles[i].Nodes[0]->GlNum;
      i2=Triangles[i].Nodes[1]->GlNum;
      i3=Triangles[i].Nodes[2]->GlNum;

      cblas_dcopy(3,Nodes[i1].X,1,X+3*i,1);
      cblas_daxpy(3,1.0,Nodes[i2].X,1,X+3*i,1);
      cblas_daxpy(3,1.0,Nodes[i3].X,1,X+3*i,1);
      alpha=1.0/3.0;
      cblas_dscal(3,alpha,X+3*i,1);

      G[i]=Triangles[i].Area;

    } 
/* 
  Set ACA Parameters
*/
  ParamACA.NCluMin=NCluMin;
  ParamACA.MaxRank=MaxRank;
  ParamACA.ACAEta=ACAEta;
  ParamACA.ACAAlpha=ACAAlpha;
  ParamACA.ACABeta=ACABeta;
  ParamACA.ACAGamma=ACAGamma;
  ParamACA.ACAEps=ACAEps;
  ParamACA.Symmetry=Symmetry;
  ParamACA.Info=1;
/* 
  Construct cluster tree for nodes
*/
    printf("Construct cluster tree ...\n"); 
    fflush(stdout);

    NClusters=ClusterTree(NumElements,G,X,&Clusters,&Permu,ParamACA);

    printf("  Done with %ld clusters\n\n",NClusters);
    fflush(stdout);
/* 
   Construct cluster symmetric pairs nodes-nodes
*/
    printf("Construct cluster pairs element-element ...\n"); 
    fflush(stdout);

    NPairs=SymmClusterPairs(Clusters,&Pairs,&NMax,&MMax,ParamACA);

    printf("  Done with %ld cluster pairs\n",NPairs);
    printf("  Maximal dimension NMax = %ld, MMax = %ld\n\n",NMax,MMax);
    fflush(stdout);

    ParamACA.NMax=NMax;
    ParamACA.MMax=MMax;

    XGlobal=X;

/*     printf("Generate H-Matrix ...\n");  */
/*     fflush(stdout); */

/*     GenHMatrix(NumElements,NumElements,0,0,Permu,Permu, */
/*                Clusters,Clusters, */
/*                NPairs,Pairs,&HMat, */
/*                GenElement,&ParamACA); */
/*     printf("  Done with %8.2f MB instead of regular %8.2f MB (%6.2f %%)\n",ParamACA.PCMemory,ParamACA.Memory,ParamACA.PCMemory/ParamACA.Memory*100.0); */
/*     printf("  Post compression with %8.2f MB instead of %8.2f MB (%6.2f %%)\n",ParamACA.PCMemory,ParamACA.AppMemory,ParamACA.PCMemory/ParamACA.AppMemory*100.0); */
/*     fflush(stdout); */
/*
  Get BlockHMatrix
*/
  BHMat.NRow=2*NumElements+1;
  BHMat.NColumn=2*NumElements+1;

  BHMat.NBlockRow=3;
  BHMat.NBlockColumn=3;

  BHMat.NDim=(long *)malloc(BHMat.NBlockRow*sizeof(long));
  assert(BHMat.NDim != NULL);

  BHMat.MDim=(long *)malloc(BHMat.NBlockColumn*sizeof(long));
  assert(BHMat.MDim != NULL);
 
  BHMat.NDim[0]=0;
  BHMat.NDim[1]=NumElements;
  BHMat.NDim[2]=2*NumElements;

  BHMat.MDim[0]=0;
  BHMat.MDim[1]=NumElements;
  BHMat.MDim[2]=2*NumElements;

  BHMat.HBlocks=(HMatrix *)malloc(BHMat.NBlockRow*BHMat.NBlockColumn*sizeof(HMatrix));
  assert(BHMat.HBlocks != NULL);
  /* BHMat.NRow=NumElements; */
  /* BHMat.NColumn=NumElements; */

  /* BHMat.NBlockRow=1; */
  /* BHMat.NBlockColumn=1; */

  /* BHMat.NDim=(long *)malloc(BHMat.NBlockRow*sizeof(long)); */
  /* assert(BHMat.NDim != NULL); */

  /* BHMat.MDim=(long *)malloc(BHMat.NBlockColumn*sizeof(long)); */
  /* assert(BHMat.MDim != NULL); */
 
  /* BHMat.NDim[0]=0; */

  /* BHMat.MDim[0]=0; */

  /* BHMat.HBlocks=(HMatrix *)malloc(BHMat.NBlockRow*BHMat.NBlockColumn*sizeof(HMatrix)); */
  /* assert(BHMat.HBlocks != NULL); */
/*
  Block 1,1
*/
    printf("Generate H-Matrix ...\n");
    fflush(stdout);

    GenHMatrix(NumElements,NumElements,0,0,Permu,Permu,
               Clusters,Clusters,
               NPairs,Pairs,BHMat.HBlocks+0,
               GenElement,&ParamACA);
    printf("  Done with %8.2f MB instead of regular %8.2f MB (%6.2f %%)\n",ParamACA.PCMemory,ParamACA.Memory,ParamACA.PCMemory/ParamACA.Memory*100.0);
    printf("  Post compression with %8.2f MB instead of %8.2f MB (%6.2f %%)\n",ParamACA.PCMemory,ParamACA.AppMemory,ParamACA.PCMemory/ParamACA.AppMemory*100.0);
    fflush(stdout);
/*
  Block 2,1
*/
    BHMat.HBlocks[1].NRow=NumElements;
    BHMat.HBlocks[1].NColumn=NumElements;
    BHMat.HBlocks[1].PermuRow=Permu;
    BHMat.HBlocks[1].PermuColumn=Permu;
    BHMat.HBlocks[1].Symmetry=ParamACA.Symmetry;
    BHMat.HBlocks[1].NBlocks=1;
    BHMat.HBlocks[1].Precond=None;
    BHMat.HBlocks[1].CBlocks=(CBlock *)malloc(1*sizeof(CBlock));
    assert(BHMat.HBlocks[1].CBlocks != NULL);

    (*BHMat.HBlocks[1].CBlocks).IBlock=0;
    (*BHMat.HBlocks[1].CBlocks).JBlock=0;
    (*BHMat.HBlocks[1].CBlocks).NBlock=NumElements;
    (*BHMat.HBlocks[1].CBlocks).MBlock=NumElements;
    (*BHMat.HBlocks[1].CBlocks).Type=Null;
    (*BHMat.HBlocks[1].CBlocks).A11=NULL;
    (*BHMat.HBlocks[1].CBlocks).A21=NULL;
    (*BHMat.HBlocks[1].CBlocks).A12=NULL;
    (*BHMat.HBlocks[1].CBlocks).A22=NULL;
    (*BHMat.HBlocks[1].CBlocks).Data=NULL;
    (*BHMat.HBlocks[1].CBlocks).Rank=0;
    (*BHMat.HBlocks[1].CBlocks).UData=NULL;
    (*BHMat.HBlocks[1].CBlocks).VData=NULL;
/*
  Block 3,1
*/
    BHMat.HBlocks[2].NRow=1;
    BHMat.HBlocks[2].NColumn=NumElements;
    BHMat.HBlocks[2].PermuRow=Permu;
    BHMat.HBlocks[2].PermuColumn=Permu;
    BHMat.HBlocks[2].Symmetry=ParamACA.Symmetry;
    BHMat.HBlocks[2].NBlocks=1;
    BHMat.HBlocks[2].Precond=None;
    BHMat.HBlocks[2].CBlocks=(CBlock *)malloc(1*sizeof(CBlock));
    assert(BHMat.HBlocks[2].CBlocks != NULL);

    (*BHMat.HBlocks[2].CBlocks).IBlock=0;
    (*BHMat.HBlocks[2].CBlocks).JBlock=0;
    (*BHMat.HBlocks[2].CBlocks).NBlock=1;
    (*BHMat.HBlocks[2].CBlocks).MBlock=NumElements;
    (*BHMat.HBlocks[2].CBlocks).Type=Null;
    (*BHMat.HBlocks[2].CBlocks).A11=NULL;
    (*BHMat.HBlocks[2].CBlocks).A21=NULL;
    (*BHMat.HBlocks[2].CBlocks).A12=NULL;
    (*BHMat.HBlocks[2].CBlocks).A22=NULL;
    (*BHMat.HBlocks[2].CBlocks).Data=NULL;
    (*BHMat.HBlocks[2].CBlocks).Rank=0;
    (*BHMat.HBlocks[2].CBlocks).UData=NULL;
    (*BHMat.HBlocks[2].CBlocks).VData=NULL;
/*
  Block 1,2
*/
    BHMat.HBlocks[3].NRow=NumElements;
    BHMat.HBlocks[3].NColumn=NumElements;
    BHMat.HBlocks[3].PermuRow=Permu;
    BHMat.HBlocks[3].PermuColumn=Permu;
    BHMat.HBlocks[3].Symmetry=ParamACA.Symmetry;
    BHMat.HBlocks[3].NBlocks=1;
    BHMat.HBlocks[3].Precond=None;
    BHMat.HBlocks[3].CBlocks=(CBlock *)malloc(1*sizeof(CBlock));
    assert(BHMat.HBlocks[3].CBlocks != NULL);

    (*BHMat.HBlocks[3].CBlocks).IBlock=0;
    (*BHMat.HBlocks[3].CBlocks).JBlock=0;
    (*BHMat.HBlocks[3].CBlocks).NBlock=NumElements;
    (*BHMat.HBlocks[3].CBlocks).MBlock=NumElements;
    (*BHMat.HBlocks[3].CBlocks).Type=Null;
    (*BHMat.HBlocks[3].CBlocks).A11=NULL;
    (*BHMat.HBlocks[3].CBlocks).A21=NULL;
    (*BHMat.HBlocks[3].CBlocks).A12=NULL;
    (*BHMat.HBlocks[3].CBlocks).A22=NULL;
    (*BHMat.HBlocks[3].CBlocks).Data=NULL;
    (*BHMat.HBlocks[3].CBlocks).Rank=0;
    (*BHMat.HBlocks[3].CBlocks).UData=NULL;
    (*BHMat.HBlocks[3].CBlocks).VData=NULL;
/*
  Block 2,2
*/
    printf("Generate H-Matrix ...\n");
    fflush(stdout);

    GenHMatrix(NumElements,NumElements,0,0,Permu,Permu,
               Clusters,Clusters,
               NPairs,Pairs,BHMat.HBlocks+4,
               GenElement,&ParamACA);
    printf("  Done with %8.2f MB instead of regular %8.2f MB (%6.2f %%)\n",ParamACA.PCMemory,ParamACA.Memory,ParamACA.PCMemory/ParamACA.Memory*100.0);
    printf("  Post compression with %8.2f MB instead of %8.2f MB (%6.2f %%)\n",ParamACA.PCMemory,ParamACA.AppMemory,ParamACA.PCMemory/ParamACA.AppMemory*100.0);
    fflush(stdout);
/*
  Block 3,2
*/
    BHMat.HBlocks[5].NRow=1;
    BHMat.HBlocks[5].NColumn=NumElements;
    BHMat.HBlocks[5].PermuRow=Permu;
    BHMat.HBlocks[5].PermuColumn=Permu;
    BHMat.HBlocks[5].Symmetry=ParamACA.Symmetry;
    BHMat.HBlocks[5].NBlocks=1;
    BHMat.HBlocks[5].Precond=None;
    BHMat.HBlocks[5].CBlocks=(CBlock *)malloc(1*sizeof(CBlock));
    assert(BHMat.HBlocks[5].CBlocks != NULL);

    (*BHMat.HBlocks[5].CBlocks).IBlock=0;
    (*BHMat.HBlocks[5].CBlocks).JBlock=0;
    (*BHMat.HBlocks[5].CBlocks).NBlock=1;
    (*BHMat.HBlocks[5].CBlocks).MBlock=NumElements;
    (*BHMat.HBlocks[5].CBlocks).Type=Null;
    (*BHMat.HBlocks[5].CBlocks).A11=NULL;
    (*BHMat.HBlocks[5].CBlocks).A21=NULL;
    (*BHMat.HBlocks[5].CBlocks).A12=NULL;
    (*BHMat.HBlocks[5].CBlocks).A22=NULL;
    (*BHMat.HBlocks[5].CBlocks).Data=NULL;
    (*BHMat.HBlocks[5].CBlocks).Rank=0;
    (*BHMat.HBlocks[5].CBlocks).UData=NULL;
    (*BHMat.HBlocks[5].CBlocks).VData=NULL;
/*
  Block 1,3
*/
    BHMat.HBlocks[6].NRow=NumElements;
    BHMat.HBlocks[6].NColumn=1;
    BHMat.HBlocks[6].PermuRow=Permu;
    BHMat.HBlocks[6].PermuColumn=Permu;
    BHMat.HBlocks[6].Symmetry=ParamACA.Symmetry;
    BHMat.HBlocks[6].NBlocks=1;
    BHMat.HBlocks[6].Precond=None;
    BHMat.HBlocks[6].CBlocks=(CBlock *)malloc(1*sizeof(CBlock));
    assert(BHMat.HBlocks[6].CBlocks != NULL);

    (*BHMat.HBlocks[6].CBlocks).IBlock=0;
    (*BHMat.HBlocks[6].CBlocks).JBlock=0;
    (*BHMat.HBlocks[6].CBlocks).NBlock=NumElements;
    (*BHMat.HBlocks[6].CBlocks).MBlock=1;
    (*BHMat.HBlocks[6].CBlocks).Type=Null;
    (*BHMat.HBlocks[6].CBlocks).A11=NULL;
    (*BHMat.HBlocks[6].CBlocks).A21=NULL;
    (*BHMat.HBlocks[6].CBlocks).A12=NULL;
    (*BHMat.HBlocks[6].CBlocks).A22=NULL;
    (*BHMat.HBlocks[6].CBlocks).Data=NULL;
    (*BHMat.HBlocks[6].CBlocks).Rank=0;
    (*BHMat.HBlocks[6].CBlocks).UData=NULL;
    (*BHMat.HBlocks[6].CBlocks).VData=NULL;
/*
  Block 2,3
*/
    BHMat.HBlocks[7].NRow=NumElements;
    BHMat.HBlocks[7].NColumn=1;
    BHMat.HBlocks[7].PermuRow=Permu;
    BHMat.HBlocks[7].PermuColumn=Permu;
    BHMat.HBlocks[7].Symmetry=ParamACA.Symmetry;
    BHMat.HBlocks[7].NBlocks=1;
    BHMat.HBlocks[7].Precond=None;
    BHMat.HBlocks[7].CBlocks=(CBlock *)malloc(1*sizeof(CBlock));
    assert(BHMat.HBlocks[7].CBlocks != NULL);

    (*BHMat.HBlocks[7].CBlocks).IBlock=0;
    (*BHMat.HBlocks[7].CBlocks).JBlock=0;
    (*BHMat.HBlocks[7].CBlocks).NBlock=NumElements;
    (*BHMat.HBlocks[7].CBlocks).MBlock=1;
    (*BHMat.HBlocks[7].CBlocks).Type=Null;
    (*BHMat.HBlocks[7].CBlocks).A11=NULL;
    (*BHMat.HBlocks[7].CBlocks).A21=NULL;
    (*BHMat.HBlocks[7].CBlocks).A12=NULL;
    (*BHMat.HBlocks[7].CBlocks).A22=NULL;
    (*BHMat.HBlocks[7].CBlocks).Data=NULL;
    (*BHMat.HBlocks[7].CBlocks).Rank=0;
    (*BHMat.HBlocks[7].CBlocks).UData=NULL;
    (*BHMat.HBlocks[7].CBlocks).VData=NULL;
/*
  Block 3,3
*/
    long P1[1];
    P1[0]=0;

    BHMat.HBlocks[8].NRow=1;
    BHMat.HBlocks[8].NColumn=1;
    BHMat.HBlocks[8].PermuRow=P1;
    BHMat.HBlocks[8].PermuColumn=P1;
    BHMat.HBlocks[8].Symmetry=ParamACA.Symmetry;
    BHMat.HBlocks[8].NBlocks=1;
    BHMat.HBlocks[8].Precond=None;
    BHMat.HBlocks[8].CBlocks=(CBlock *)malloc(1*sizeof(CBlock));
    assert(BHMat.HBlocks[8].CBlocks != NULL);

    (*BHMat.HBlocks[8].CBlocks).IBlock=0;
    (*BHMat.HBlocks[8].CBlocks).JBlock=0;
    (*BHMat.HBlocks[8].CBlocks).NBlock=1;
    (*BHMat.HBlocks[8].CBlocks).MBlock=1;
    (*BHMat.HBlocks[8].CBlocks).Type=Dense;
    (*BHMat.HBlocks[8].CBlocks).A11=NULL;
    (*BHMat.HBlocks[8].CBlocks).A21=NULL;
    (*BHMat.HBlocks[8].CBlocks).A12=NULL;
    (*BHMat.HBlocks[8].CBlocks).A22=NULL;
    (*BHMat.HBlocks[8].CBlocks).Data=(double complex *)malloc(1*sizeof(double complex));
    assert((*BHMat.HBlocks[8].CBlocks).Data != NULL);
    (*BHMat.HBlocks[8].CBlocks).Data[0]=1.0+I*0.0;
    (*BHMat.HBlocks[8].CBlocks).Rank=0;
    (*BHMat.HBlocks[8].CBlocks).UData=NULL;
    (*BHMat.HBlocks[8].CBlocks).VData=NULL;
/* 
   Construct right hand side
*/
    printf("Generate RHS ...\n"); 
    fflush(stdout);

    /* xv=(double complex *)malloc(NumElements*sizeof(double complex)); */
    xv=(double complex *)malloc((2*NumElements+1)*sizeof(double complex));
    assert(xv != NULL);

    /* yv=(double complex *)malloc(NumElements*sizeof(double complex)); */
    yv=(double complex *)malloc((2*NumElements+1)*sizeof(double complex));
    assert(yv != NULL);

    /* for (i=0; i < NumElements; i++) */
    for (i=0; i < 2*NumElements+1; i++)
      xv[i]=1.0/(1.0+((double) i)*((double) i))+I*1.0/(2.0+((double) i)*((double) i));

    BlockHMatrixMult(BHMat,xv,yv);

/*     HMatrixMult(HMat,xv,yv); */

    printf("  Done ...\n"); 
    fflush(stdout);

    printf("Start Block LU Solver ...\n"); 
    fflush(stdout);

//    BlockHLUSolveMM(1,&BHMat,xv,yv,ACAEps);
///*     HLUSolveMM(1,&HMat,xv,yv,ACAEps); */
GMResEps=1.0e-4;
  GMResMaxIter=20000;

  zset(NumElements,0.0+I*0.0,x,1);

  Result=zfgmres1(NumElements,(void (*)(double complex *,double complex *,void *))BlockHMatrixMultZ,&BlockHMat,y,x,GMResEps,GMResMaxIter,&GMResMaxIter,GMResPrint);

  printf("  \n");
  printf("  Done with %d iterations\n\n",GMResMaxIter);
  fflush(stdout); 
    printf("  Done ...\n"); 
    fflush(stdout);

    for (i=0; i < 20; i++)
      printf("i,x[i]= %ld %25.15e %25.15e\n",i,creal(xv[i]),cimag(xv[i]));

      printf("\n");

/*     for (i=NumElements-20; i < NumElements; i++) */
    for (i=NumElements+1-20; i < NumElements; i++)
      printf("i,x[i]= %ld %25.15e %25.15e\n",i,creal(xv[i]),cimag(xv[i]));


/*     long ipos,jpos; */
/*     double complex aij,fronorm, cerror; */
 
/*     cerror=0.0+0.0*I; */
/*     fronorm=0.0+0.0*I; */

/*     for (jpos=0; jpos < NumElements; jpos++) */
/*       { */
/*         zset(NumElements,0.0+0.0*I,xv,1); */
/* 	xv[jpos]=1.0+0.0*I; */
/*         HMatrixMult(HMat,xv,yv); */
/* 	for (ipos=0; ipos < NumElements; ipos++) */
/* 	  { */
/* 	    aij=GenElement(ipos,jpos); */
/*             fronorm+=aij*conj(aij); */
/* 	    cerror+=(yv[ipos]-aij)*conj(yv[ipos]-aij); */
/* 	  } */
/*       } */

/*     printf("Fro-Norm = %25.15e Rel-Error = %25.15e\n",sqrt(creal(fronorm)),sqrt(creal(cerror/fronorm))); */

/*
  Single elements
*/
/*     while(1 == 1) */
/*       { */
/* 	printf("i="); */
/* 	scanf("%ld",&ipos); */
/*         if (ipos == -1) break; */
/* 	printf("j="); */
/* 	scanf("%ld",&jpos); */
/* 	printf("i,j= %ld %ld\n",ipos,jpos); */
        
/*         zset(NumElements,0.0+0.0*I,xv,1); */
/* 	xv[jpos]=1.0+0.0*I; */

/*         HMatrixMult(HMat,xv,yv); */

/*         printf("App = %14.8e + %14.8e*I\nExa = %14.8e + %14.8e*I\n",creal(yv[ipos]),cimag(yv[ipos]),creal(GenElement(ipos,jpos)),cimag(GenElement(ipos,jpos))); */
/*       } */
/*
  Free memory
*/
/*     free(Nodes[0].Nb); free(Nodes[0].KNb); */
    free(Nodes); free(Triangles);
    free(G); free(X); free(Clusters); free(Pairs);
/*     FreeHMatrix(&HMat); */
    /* FreeBlockHMatrix(&BHMat); */
    free(xv); free(yv);
}
