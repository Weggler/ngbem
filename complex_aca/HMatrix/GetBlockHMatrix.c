//
//  routine name  -  get_block_hmatrix
//
//------------------------------------------------------------------
//
//  last revision -  Sep 12
//
//  purpose       -  generate block hierarchic matrix
//
//  in                                                              
//        NRBLOCK - number of block
//      NROW,NCOL - row, column number
//       NKIND1,2 - pairing number and type clusters
//            0,0 == medg-medg
//            1,1 == mdle-mdle
//            0,1 == medg-mdle
//            1,0 == mdle-medg                                      
//     NDIM,MDIM  - block dimension
//     IPOS,JPOS  - position of first element of the block matrix
//                  in the global matrix, starting with 0
//
//  out          - generated block matrix in (*BHMAT)
//
// GLOBAL VARIABLES: *BHMAT, CLUSTERS, NPAIRS, PAIRS, PARAMACA 
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/GlobalVariables.h"

double complex GenElement(long IEl,long JEl);

#define get_block_hmatrix get_block_hmatrix_

void get_block_hmatrix(int *NRBLOCK, 
                       int *NROW,   int *NCOL,
                       int *NKIND1, int *NKIND2,
                       int *NDIM,   int *MDIM, 
                       int *IPOS,   int *JPOS)
{
   char temp[10]; 
   char file[10]; 
   int iprint = 1;
   long nrblock;
   long np;
   long iclu,jclu;
   long nrow,ncol;
   long ndim,mdim;
   long ipos,jpos;
   int FileD;

   nrblock = (long) *NRBLOCK;
   nrow = (long) *NROW;
   ncol = (long) *NCOL;
   iclu = (long) *NKIND1;
   jclu = (long) *NKIND2;
   ndim = (long) *NDIM;
   mdim = (long) *MDIM;
   ipos = (long) *IPOS;
   jpos = (long) *JPOS;

   if(iclu == 0 && jclu == 0)      // medg-medg
   {
     np = 0;
   }
   else if(iclu == 1 && jclu == 1) // mdle-mdle
   {
     np = 1;
   }
   else if(iclu == 0 && jclu == 1) // medg-mdle
   {
     np = 2;
   }
   else if(iclu == 1 && jclu == 0) // mdle-medg
   {
     np = 3;
   }

   /* Parameters */
   PARAMACA.NMax=NMAX[iclu];
   PARAMACA.MMax=MMAX[jclu];
   
   /* Block position in (*BHMAT) */
   (*BHMAT).NDim[nrow] = ipos;
   (*BHMAT).MDim[ncol] = jpos;

   sprintf( temp, "mat%ld.bin", nrblock ); 
   strcpy( file, temp ); 
   if (access(file, F_OK) == 0)
   {
     switch(PARAMACA.Symmetry)
     {
       case NonSymmetric:
         LoadHMatrix(((*BHMAT).HBlocks+nrblock),file);
         break;
       case ComplexSymmetric:
         LoadCSymHMatrix(((*BHMAT).HBlocks+nrblock),file);
         break;
     }
   }  
   else 
   {
     GenHMatrix(ndim,mdim,ipos,jpos,
                PERMU[iclu],PERMU[jclu],
                CLUSTERS[iclu],CLUSTERS[jclu],
                NPAIRS[np],PAIRS[np],(*BHMAT).HBlocks+nrblock,
                GenElement,&PARAMACA);
     if(iprint == 1) 
     {
       printf("GetBlockHMatrix\n");
       printf("  Matrix Block %ld\n", nrblock);
       printf("  Position     %ld, %ld\n", ipos, jpos);
       printf("  Done with %8.2f MB instead of regular %8.2f MB (%6.2f %%)\n",
              PARAMACA.PCMemory,PARAMACA.Memory,
              PARAMACA.PCMemory/PARAMACA.Memory*100.0);
       printf("  Post compression with %8.2f MB instead of %8.2f MB (%6.2f %%)\n",
              PARAMACA.PCMemory,PARAMACA.AppMemory,
              PARAMACA.PCMemory/PARAMACA.AppMemory*100.0);
     }
     switch(PARAMACA.Symmetry)
     {
       case NonSymmetric:
         SaveHMatrix(*((*BHMAT).HBlocks+nrblock),file);
         break;
       case ComplexSymmetric:
         SaveCSymHMatrix(*((*BHMAT).HBlocks+nrblock),file);
         break;
     }
/* NOTE: to run par.sh you should be sure that writing all mat%ld.bin files is  
         fully completed before calling maxwell_parallel.
         this is why a dummy file mat%ld.done is generated per HBlock 
         (their total number is read by par.sh) */
     sprintf(temp, "mat%ld.done", nrblock ); 
     FileD = open(temp,O_WRONLY | O_CREAT, S_IRWXU); assert(FileD != -1); close(FileD);

/* NOTE: enable the exit command when compiling the block generation for 
         to run aca in 'parallel mode' otherwise comment it out */
//	 exit(0);
   }
}
