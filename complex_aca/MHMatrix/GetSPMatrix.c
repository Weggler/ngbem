//------------------------------------------------------------------
//                                                                  
//  routine name   -  cG_to_dcG_interface
//                                                                  
//------------------------------------------------------------------
//
//  last revision  -  Feb 13
//
//  purpose        -  generate global connectivity data called 
//                    by multiplication routine in GMRes solver
//
//  in:
//       Ndoftot   -  number of cG dof's
//       Ndofmaca  -  number of gcG dof's
//       Nr_con    -  list of number of gcG dofs per cG dof
//       Nodes_con -  corresponding list of nodes in MACA numbering 
//
//  out:          
//
// GLOBAL VARIABLES: *NR_CON, *NODES_CON, NDOFTOT, NDOFMACA
//
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/GlobalVariables.h"

#define cg_to_dcg_interface cg_to_dcg_interface_

void cg_to_dcg_interface(int* Ndoftot,int* Ndofmaca, int* Nr_con, int* Nodes_con)
{
   /* local variables */
   int k;
   int l;
   int N = 0;
   NDOFTOT = *Ndoftot;
   NDOFMACA = *Ndofmaca;

   NR_CON = (int*) malloc(NDOFTOT*sizeof(int));
   for(k=0; k<NDOFTOT; k++)
   {
     NR_CON[k] = Nr_con[k];
     N += Nr_con[k];
     
   }
   assert(N == *Ndofmaca);
   NODES_CON = (int*) malloc(N*sizeof(int));
   N = 0;
   for(k=0; k<NDOFTOT; k++)
   {
     for(l=0; l<NR_CON[k]; l++)
     {
       NODES_CON[N] = Nodes_con[k+l*NDOFTOT]-1;
       N++;
     }
   }
}
