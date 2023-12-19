/*******************************************************************\
*                                                                   *
*                          S. Rjasanow C-Software                   *
*                                                                   *
*                       Adaptive Cross Approximation                *
*                                                                   *
*                                20.10.2010                         *
\*******************************************************************/
/*
  Last change 28.10.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void aca_connected_dof_(int *, int *, int *, int *);
void aca_dof_(int *, int *);

void BlockPrePrecond(BlockHMatrix *BlockHMat,double complex (*ElementName)(long,long),double complex *y)
{
/*
  Local variables
*/
  long i,j,k,l;
  long ib,jb,N,M,NB,MB,NBlock,MBlock,IPos,JPos,Error,NumNB,NumkNB;
  long kb,ke,im,jm,IEl,JEl;
  long *ListNB,*ListkNB;
  PrecondType Precond;
  double complex *XCopy,*YCopy,*Mat,*xn,*yn;
  HMatrix Block;
  SparseRow *SparseRows;
/*
  Feko variables
*/
  int iF,numF;
  int *listF;
/*
  LAPACK variables
*/
//  MKL_INT NL,NRHSL,LDAL,LDBL,INFOL;
//  MKL_INT *IPIVL;
//  MKL_Complex16 *MATL,*YL;
  int NL,NRHSL,LDAL,LDBL,INFOL;
  int *IPIVL;
  double complex *MATL,*YL;

  N=(*BlockHMat).NRow;
  M=(*BlockHMat).NColumn;
  
  NB=(*BlockHMat).NBlockRow;
  MB=(*BlockHMat).NBlockColumn;

  Precond = (*BlockHMat).Precond;

  switch(Precond)
  {
    case None:
      break;

    case Jacobi:
  
      BlockHMat->PrecondIData=NULL;
      BlockHMat->PrecondDData=(double complex *)malloc(N*sizeof(double complex));
  
      for (i=0; i < N; i++)
      {
        (BlockHMat->PrecondDData)[i]=ElementName(i,i);
        y[i]=y[i]/(BlockHMat->PrecondDData)[i];
      }
      break;
    
	case TriDiag:
  
      BlockHMat->PrecondIData=NULL;
      BlockHMat->PrecondDData=(double complex *)malloc(3*N*sizeof(double complex));
      assert(BlockHMat->PrecondDData != NULL);
  
      (BlockHMat->PrecondDData)[0]=0.0+0.0*I;
      (BlockHMat->PrecondDData)[N+0]=ElementName(0,1);
      (BlockHMat->PrecondDData)[2*N+0]=ElementName(0,0);
  
      for (i=1; i < N-1; i++)
      {
        (BlockHMat->PrecondDData)[i]=ElementName(i-1,i);
        (BlockHMat->PrecondDData)[N+i]=ElementName(i,i+1);
        (BlockHMat->PrecondDData)[2*N+i]=ElementName(i,i);
      }
  
      (BlockHMat->PrecondDData)[N-1]=ElementName(N-2,N-1);
      (BlockHMat->PrecondDData)[2*N-1]=0.0+0.0*I;
      (BlockHMat->PrecondDData)[3*N-1]=ElementName(N-1,N-1);
  
      SolveTriDiag(N-1,(BlockHMat->PrecondDData),(BlockHMat->PrecondDData)+N,(BlockHMat->PrecondDData)+2*N,y,y);       
      break;
    
    case ILU:
    
//      printf("In ILU:\n");
//    
//      BlockHMat->PrecondIData=NULL;
//      BlockHMat->PrecondDData=NULL;
//    
//      GetBlockSparseStructure(BlockHMat);
//    
//      XCopy=(double complex *)malloc(N*sizeof(double complex));
//      assert(XCopy != NULL);
//    
//      YCopy=(double complex *)malloc(N*sizeof(double complex));
//      assert(YCopy != NULL);
///*
//  Permutation of the given vector with P_Row^T
//*/
//      for (ib=0; ib < NB; ib++)
//      {
//        IPos=(*BlockHMat).NDim[ib];
//    
//        Block=(*BlockHMat).HBlocks[ib];
//        NBlock=Block.NRow;
//    
//        for (i=0; i < NBlock; i++)
//        {
//          YCopy[IPos+i]=y[IPos+Block.PermuRow[i]];
//        }
//      }
///*
//  Incomplete LU
//*/
//      Error=zlusolC(YCopy,XCopy,BlockHMat->lu);
//      assert(Error == 0);
///*
//  Permutation of the resulting vector with P_Column^T
//*/
//      for (jb=0; jb < NB; jb++)
//      {
//        JPos=(*BlockHMat).MDim[jb];
//      
//        Block=(*BlockHMat).HBlocks[jb*NB];
//        MBlock=Block.NColumn;
//      
//        for (j=0; j < MBlock; j++)
//          y[JPos+Block.PermuColumn[j]]=XCopy[JPos+j];
//      }
///*
//  Free memory
//*/          
//      free(XCopy); free(YCopy);
//      break;

    case BasisElement:

      printf("Start Preconditioning ...\n"); 
      fflush(stdout);

      BlockHMat->PrecondIData=NULL;
      BlockHMat->PrecondDData=NULL;
      
      SparseRows=(SparseRow *)malloc(N*sizeof(SparseRow));
      assert(SparseRows != NULL);
      
      listF=(int *)malloc(N*sizeof(int));
      assert(listF != NULL);
      
      ListkNB=(long *)malloc(N*sizeof(long));
      assert(ListkNB != NULL);

// ...start with local elements all with global numbers starting with 0
      for(i=0; i < N; i++)
      {
         NumNB=0;
         JoinSets(1,&i,&NumNB,&ListNB);
         ke=0;
         for(j=0; j <= LevelNB; j++)
         {
           kb=ke;
           ke=NumNB;
           for(k=kb; k < ke; k++)
    	   {
// ..........change C-numbering to Fortran-numbering
    	     iF=(int)ListNB[k]+1;
             aca_connected_dof_(&iF,(int*) &N, &numF,listF);

             NumkNB=(long)numF;
// ..........change C-numbering to Fortran-numbering
             for (l=0; l < NumkNB; l++)
               ListkNB[l]=(long)listF[l]-1; 
             JoinSets(NumkNB,ListkNB,&NumNB,&ListNB);
    	   }
         }

// ......concatenate ListNB with basis elements (*BlockHMat).BE 
//       leading to new ListNB of NumNB elements. 
//       ListNB == list of all dofs that are used to set up the 
//                 preconditioning matrix, 
//                 the dof-numbering is global starting with 0
//        JoinSets((*BlockHMat).NumBE,(*BlockHMat).BE,&NumNB,&ListNB);

        Mat=(double complex *)malloc(NumNB*NumNB*sizeof(double complex));
        assert(Mat != NULL);
        
        yn=(double complex *)malloc(NumNB*sizeof(double complex));
        assert(yn != NULL);

        IPos=0;
        for(jm=0; jm < NumNB; jm++)
        {                     
          for(im=0; im < NumNB; im++)
          {
            JEl=ListNB[jm];     
            IEl=ListNB[im]; // in C-numbering
            aca_dof_((int*) &IEl,(int*) &JEl);
            Mat[IPos]=GenElement(IEl,JEl);
            IPos++;
          }
        }

        for(im=0; im < NumNB; im++)
          if (ListNB[im] == i)
            yn[im]=1.0+I*0.0;
          else
            yn[im]=0.0+I*0.0;
/*
  LAPACK call
*/
//          IPIVL=(MKL_INT *)malloc(NumNB*sizeof(MKL_INT));

//          NL=(MKL_INT)NumNB;
//          NRHSL=1;
//          LDAL=(MKL_INT)NumNB;
//          LDBL=(MKL_INT)NumNB;
// 
//          MATL=(MKL_Complex16 *)Mat;
//          YL=(MKL_Complex16 *)yn;

        IPIVL=(int *)malloc(NumNB*sizeof(int));
        assert(IPIVL != NULL);
        
        NL=(int)NumNB;
        NRHSL=1;
        LDAL=(int)NumNB;
        LDBL=(int)NumNB;
        
        MATL=(double complex *)Mat;
        YL=(double complex *)yn;
        
        zgesv_(&NL,&NRHSL,MATL,&LDAL,IPIVL,YL,&LDBL,&INFOL);
        if(INFOL != 0)
        {
          printf("BlockPrePrecond: %ld\n\n", i);
          fflush(stdout);
        }
        assert(INFOL == 0);
        
        SparseRows[i].NI=NumNB;
        SparseRows[i].J=ListNB;
        SparseRows[i].AIJ=yn;
/* 
  Free memory
*/         
        free(MATL); free(IPIVL);

        printf("   %6.2f %% done\r",((double)i)/((double)N)*100.0);
        fflush(stdout);
     }

     BlockHMat->PrecondData=SparseRows;
     
     XCopy=(double complex *)malloc(N*sizeof(double complex));
     assert(XCopy != NULL);
     
     cblas_zcopy(N,y,1,XCopy,1);
     
     SparseRowMult(N,SparseRows,XCopy,y);
/*
  Free memory
*/          
     free(XCopy); free(listF); free(ListkNB);
    break;
  }

  printf("\n");
  fflush(stdout);
}

/*****************************************************************/
/* 
   Test the preconditioning for higher order meshes:
   IDEA: if LevelNB is high enough, then for each i
         in the loop the full matrix is generated and the system
         is solved by the direct solver zgesv_
      => turn off the gmres solver just analyse what you get from the 
         "predended" preconditioning which is in fact solving now... as
         you should obtain the i-th component of the soution
         vector which is to double check with the direct solver 
         aca version. 
   IDEA: to exclude errors in aca_connected_dof_ I simply substituted
         the double loop by the following lines, i.e., I set up a 
         simple run through all dofs ordered from 0 to NumNB-1 !
*/
//        IPos=0;
//        printf("BlockPrePrecond: %ld %ld\n", NumNB, N);
//        for(jm=0; jm < NumNB; jm++)
//        {                     
//          for(im=0; im < NumNB; im++)
//          {
//// IMPORTANT: reset JEl here because it got possibly
////            changed in last aca_basis_elements_ call ! 
//            JEl=jm;     
//            IEl=im;         // in C-numbering
//            ListNB[jm]=jm;  // Fortran -> C-numbering 
//            aca_dof_((int*) &IEl,(int*) &JEl);
//            Mat[IPos]=GenElement(IEl,JEl);
//            IPos++;
//          }
//        }

/*****************************************************************/

/* first I thought that we should also add the neighbors of the 
   basis elements... if so, then insert the following lines after calling
  "JoinSets((*BlockHMat).NumBE,(*BlockHMat).BE,&NumNB,&ListNB);"
*/

//        ke  = NumNB;
//        for (l=0; l <ke; l++)
//        {
//          printf("BlockPrePrecond before: %ld, %ld\n", l, ListNB[l]);
//          fflush(stdout);
//        }
//        for(k=0; k < ke; k++)
//        {
//    	    iF=(int)ListNB[k]+1;
//          aca_connected_dof_(&iF,(int*) &N, &numF,listF);
//              
//          NumkNB=(long)numF;
//           
//          for (l=0; l < NumkNB; l++)
//          {
//            ListkNB[l]=(long)listF[l]-1; 
//          }
//          
//          JoinSets(NumkNB,ListkNB,&NumNB,&ListNB);
//        }
//        for (l=0; l < NumNB; l++)
//        {
//          printf("BlockPrePrecond after: %ld, %ld\n", l, ListNB[l]);
//          fflush(stdout);
//        }
//

/*****************************************************************/
        
