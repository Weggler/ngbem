/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                14.01.2011                                   *
\*****************************************************************************/
/*
  Last change 14.01.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HDToA(CBlock *A,double HLUEps)
{
  HCheckM(A);
/*
  Local variables
*/
  long NA,MA;
  long i,j,Rank;
  double FroNorm,Crit,Sum;
  double complex *U,*V;
/*
  LAPACK variables
*/
//  MKL_INT NL,ML,LWorkL,Info;
//  MKL_Complex16 *AL,*UL,*VL,*WorkL;
  int NL,ML,LWorkL,Info;
  double complex *AL,*UL,*VL,*WorkL;
  double *S,*RWork;
/*
  Dimension of the matrix A
*/
  NA=(*A).NBlock;
  MA=(*A).MBlock;

  if (lmin(MA,NA) == 1) 
    return;

  assert((*A).Type == Dense);

//  NL=(MKL_INT)NA;
//  ML=(MKL_INT)MA;
//  AL=(MKL_Complex16 *)malloc(NA*MA*sizeof(MKL_Complex16));
  NL=(int)NA;
  ML=(int)MA;
  AL=(double complex *)malloc(NA*MA*sizeof(double complex));
  assert(AL != NULL);

  cblas_zcopy(NA*MA,(*A).Data,1,AL,1);

//  UL=(MKL_Complex16 *)malloc(NA*NA*sizeof(MKL_Complex16));
  UL=(double complex *)malloc(NA*NA*sizeof(double complex));
  assert(UL != NULL);

//  VL=(MKL_Complex16 *)malloc(MA*MA*sizeof(MKL_Complex16));
  VL=(double complex *)malloc(MA*MA*sizeof(double complex));
  assert(VL != NULL);

//  WorkL=(MKL_Complex16 *)malloc(10*lmax(NA,MA)*sizeof(MKL_Complex16));
  WorkL=(double complex *)malloc(10*lmax(NA,MA)*sizeof(double complex));
  assert(WorkL != NULL);

//  LWorkL=10*(MKL_INT)lmax(NA,MA);
  LWorkL=10*(int)lmax(NA,MA);

  S=(double *)malloc(lmax(NA,MA)*sizeof(double));
  assert(S != NULL);

  RWork=(double *)malloc(5*lmin(NA,MA)*sizeof(double));
  assert(RWork != NULL);

//  ZGESVD_("S","S",&NL,&ML,AL,&NL,S,UL,&NL,VL,&ML,WorkL,&LWorkL,RWork,&Info);
  zgesvd_("S","S",&NL,&ML,AL,&NL,S,UL,&NL,VL,&ML,WorkL,&LWorkL,RWork,&Info);
  assert(Info == 0);

  FroNorm=0.0;

  for (i=0; i < lmin(NA,MA); i++)
    FroNorm+=S[i]*S[i];

  Crit=HLUEps*HLUEps*FroNorm;
/*
  Compression check
*/
  Sum=0.0;
  i=lmin(NA,MA)-1;

  do
    {
      Rank=i+1;
      Sum+=S[i]*S[i];
      i-=1;
    }
  while ((Sum < Crit) && (Rank > 0));
/*
  Successful compression ?
*/
  if (Rank*(NA+MA) < NA*MA)
    {
      (*A).Type=Admissible;

      free((*A).Data);

      (*A).Data=NULL;

      (*A).Rank=Rank;

      (*A).UData=(double complex *)malloc(NA*Rank*sizeof(double complex));
      assert((*A).UData != NULL);

      (*A).VData=(double complex *)malloc(MA*Rank*sizeof(double complex));
      assert((*A).UData != NULL);

      U=(double complex *)UL;
      V=(double complex *)VL;

      for (j=0; j < Rank; j++)
      for (i=0; i < NA; i++)
	(*A).UData[j*NA+i]=S[j]*U[j*NA+i];

      for (j=0; j < Rank; j++)
      for (i=0; i < MA; i++)
	(*A).VData[j*MA+i]=conj(V[i*MA+j]);
    }
/*
  Free memory
*/
  free(AL); free(UL); free(VL); free(WorkL); free(S); free(RWork);
      
  HCheckM(A);
}
