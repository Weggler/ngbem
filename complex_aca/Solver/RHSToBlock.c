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

void RHSToBlock(long NRHS,BlockHMatrix *BHMat,double complex *Y,CBlock *BY)
{
  /* Local variables */
  long N,M,NB,MB;
  long N2,ND;
  long jr,ib,i,ia,ie;

  /* Block dimension of the matrix */
  NB=(*BHMat).NBlockRow; 
  MB=(*BHMat).NBlockColumn;
  assert(NB == MB);

  /* Dimension of the matrix */
  N=(*BHMat).NRow; 
  M=(*BHMat).NColumn;
  assert(N == M);

  /* Next power of 2 */
  N2=(long)pow(2,(long)(log(NB-0.5)/log(2.0)+1));

  ND=N2-NB;

  (*BY).IBlock=0;
  (*BY).JBlock=0;
  (*BY).NBlock=N+ND;
  (*BY).MBlock=NRHS;
  printf("RHSToBLock NB=%d, MB=%d, N=%d, M=%d, B.NBlock=%d, B.MBlock=%d\n", NB, MB, N, M, N+ND, NRHS);

  (*BY).Type=Dense;

  (*BY).A11=NULL;
  (*BY).A21=NULL;
  (*BY).A12=NULL;
  (*BY).A22=NULL;

  (*BY).Rank=0;
  (*BY).UData=NULL;
  (*BY).VData=NULL;

  (*BY).Data=(double complex *)malloc((N+ND)*NRHS*sizeof(double complex));
  assert((*BY).Data != NULL);
 
  /* Permutation of the RHS */
  for (jr=0; jr < NRHS; jr++)
  {
    for (ib=0; ib < NB; ib++)
    {
      ia=(*BHMat).NDim[ib];
      if(ib == NB-1)
        ie=N;
      else
        ie=(*BHMat).NDim[ib+1];

//      printf("RHSToBlock, ie, ia: %d %d %d\n", (*BHMat).NDim[ib+1], ie, ia);
	  for (i=ia; i < ie; i++)
      {
//         printf("RHSToBlock, Test: %d %d %d\n",i,i-ia,((*BHMat).HBlocks[ib]).PermuRow[i-ia]);
         (*BY).Data[jr*(N+ND)+i]=Y[jr*N+ia+((*BHMat).HBlocks[ib]).PermuRow[i-ia]];
      }
    }
  }
  /* Additional zeros */
  for (jr=0; jr < NRHS; jr++)
  {
    for (i=0; i < ND; i++)
    {
      (*BY).Data[jr*(N+ND)+N+i]=0.0+I*0.0;
    }
  }
}
