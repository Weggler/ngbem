//------------------------------------------------------------------
//                                                                  
//  routine name  -  get_mhmatrix
//                                                                  
//------------------------------------------------------------------
//
//  last revision - Feb 13
//
//  purpose       - generate (block) matrix-driven H-matrix
//
//  in:
//        Nrblock - number of block
//      Nrow,Ncol - number of block in block matrix
//       Nkind1,2 - pairing number and type clusters
//      Ndim,Mdim - dimension of global matrix
//      Ipos,Jpos - position of first element of the block matrix
//                  in the global matrix, starting with 0
//
//  out:          
//          BHMAT- block HMatrix  (global variable)
//
// GLOBAL VARIABLES: *BHMAT, CLUSTERS, NPAIRS, PAIRS, PARAMACA, NLOC
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/GlobalVariables.h"

void GenMElement(long IEl,long JEl,long IDim,long JDim, 
                 double complex *Zloc);

#define get_mhmatrix get_mhmatrix_

void get_mhmatrix(int *Nrblock, 
                  int *Nrow,   int *Ncol,
                  int *Nkind1, int *Nkind2,
                  int *Ndim,   int *Mdim, 
                  int *Ipos,   int *Jpos)
{
/* 
  local variables
*/
   char temp[10]; 
   char file[10]; 
   int  iprint = 1;
   long i;
   long help;
   long nrblock;
   long np;
   long iclu,jclu;
   long nrow,ncol;
   long ndim,mdim;
   long ipos,jpos;

   nrblock = (long) *Nrblock; // not used
   nrow = (long) *Nrow;       // not used
   ncol = (long) *Ncol;       // not used
   iclu = (long) *Nkind1;
   jclu = (long) *Nkind2;
   ndim = (long) *Ndim;
   mdim = (long) *Mdim;
   ipos = (long) *Ipos;
   jpos = (long) *Jpos;

   if(iclu == 0 && jclu == 0) // mdle-mdle
     np = 0;
/* 
  Maximal cluster pair dimensions 
*/
   PARAMACA.NMax=0;
   PARAMACA.MMax=0;
   for(i=0;i<CLUSTERS[iclu]->Number;i++)
   {
	 help = NLOC[iclu][PERMU[iclu][i]+1]-NLOC[iclu][PERMU[iclu][i]];
     if(help > PARAMACA.NMax)
       PARAMACA.NMax = help;
   }    
   PARAMACA.NMax*=NMAX[iclu];

   for(i=0;i<CLUSTERS[jclu]->Number;i++)
   {
     help = NLOC[jclu][PERMU[jclu][i]+1]-NLOC[jclu][PERMU[jclu][i]];
     if(help > PARAMACA.MMax)
       PARAMACA.MMax = help;
   }    
   PARAMACA.MMax*=MMAX[jclu];
/* 
  Block position in (*BHMAT) 
*/
   (*BHMAT).NDim[nrow] = ipos;
   (*BHMAT).MDim[ncol] = jpos;

   GenMHMatrix(NLOC[iclu],NLOC[jclu],
               ndim,mdim,
               ipos,jpos,
               PERMU[iclu],PERMU[jclu],
               CLUSTERS[iclu],CLUSTERS[jclu],
               NPAIRS[np],PAIRS[np],
               (*BHMAT).HBlocks+nrblock,
               GenMElement,
               &PARAMACA);
/* 
  Test prints
*/
   if(iprint == 1) 
   {
     printf("  Matrix Block  %ld\n", nrblock);
     printf("  Position      %ld, %ld\n", ipos, jpos);
//     printf("  Cluster pairs %ld x %ld\n",
//            CLUSTERS[iclu]->Number,CLUSTERS[jclu]->Number);
//     printf("  Max dim of cluster pair matrices %ld x %ld\n",
//            PARAMACA.NMax,PARAMACA.MMax);
//     printf("\n");
     printf("  Done with %8.2f MB instead of regular %8.2f MB (%6.2f %%)\n",
            PARAMACA.PCMemory,PARAMACA.Memory,
            PARAMACA.PCMemory/PARAMACA.Memory*100.0);
     printf("  Post compression with %8.2f MB instead of %8.2f MB (%6.2f %%)\n",
            PARAMACA.PCMemory,PARAMACA.AppMemory,
            PARAMACA.PCMemory/PARAMACA.AppMemory*100.0);
     fflush(stdout);
   }
}
