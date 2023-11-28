/*------------------------------------------------------------------*/
/*                                                                  */
/* S. Rjasanow: Adaptive Cross Approximation                        */
/*                                                                  */
/*  routine name  -  aca_init_precond                               */
/*                                                                  */
/*------------------------------------------------------------------*/
/*                                                                  */
/*  last revision -  Nov 11                                         */
/*  purpose       -  init preconditioning of block H-matrix         */
/*                                                                  */
/*  in                                                              */
/*        NRBLOCK - number of block                                 */
/*         NKIND1 - first part that indicates the cluster type
              0,0 == medg-medg
              1,1 == mdle-mdle
              0,1 == medg-mdle                                      */
/*                                                                  */
/*  out                                                             */
/*                                                                  */
/*  GLOBAL VARIABLES: *BHMAT, CLUSTERS, LevelBE                     */
/*                                                                  */
/*------------------------------------------------------------------*/

/* Includes for all functions */
# include "../CInclude/GlobalVariables.h"

#define aca_init_precond aca_init_precond_
void aca_init_precond(int *NRBLOCK, int *NKIND1)
{
   long nrblock;
   long np;
   long jclu,iclu;

   nrblock = (long) *NRBLOCK;
   iclu    = (long) *NKIND1;

   /* get value of precondition flag of (*BHMAT) */
   Precond = (*BHMAT).Precond;

   /* Get global basis elements */
   if(Precond == BasisElement)
   {
      for (jclu=0; jclu < NCLUSTERS[iclu]; jclu++)
      {
         if((*(CLUSTERS[iclu]+jclu)).Level == LevelBE)
         {
           JoinSets(1,&((*(CLUSTERS[iclu]+jclu)).BasisElement),&((*BHMAT).NumBE),&((*BHMAT).BE));
         }
      }
   }
}
