/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                08.05.2009                                   *
\*****************************************************************************/
/*
  Last change 08.02.2010
*/

/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void PostCompression(long N,long M,long *aRank,
                     double complex *aU,double complex *aV,
                     double Eps)
{
/*
  Local variables
*/
  long i,j,k,Rank,NewRank;
  double FroNorm,Crit,Sum,*S,*RWork;
  double complex *U,*V,*TauU,*TauV,*Work,*R,*UR,*VR;
/*
  LAPACK variables
*/
  int NL,ML,RankL,LWork,Info;
  double complex *UL,*VL,*TauUL,*TauVL,*WorkL,*RL,*URL,*VRL;
  double complex *aUL,*aVL;
//  MKL_Complex16 *UL,*VL,*TauUL,*TauVL,*WorkL,*RL,*URL,*VRL;
//  MKL_Complex16 *aUL,*aVL;
/*
  Old rank
*/
  Rank=*aRank;
/*
  Allocate memory
*/
  U=(double complex *)malloc(N*Rank*sizeof(double complex));
  assert(U != NULL);

  V=(double complex *)malloc(M*Rank*sizeof(double complex));
  assert(V != NULL);

  TauU=(double complex *)malloc(lmin(N,Rank)*sizeof(double complex));
  assert(TauU != NULL);

  TauV=(double complex *)malloc(lmin(M,Rank)*sizeof(double complex));
  assert(TauV != NULL);

  Work=(double complex *)malloc(10*lmax(N,M)*sizeof(double complex));
  assert(Work != NULL);

  RWork=(double *)malloc(5*Rank*sizeof(double));
  assert(RWork != NULL);

  R=(double complex *)malloc(Rank*Rank*sizeof(double complex));
  assert(R != NULL);

  UR=(double complex *)malloc(Rank*Rank*sizeof(double complex));
  assert(UR != NULL);

  VR=(double complex *)malloc(Rank*Rank*sizeof(double complex));
  assert(VR != NULL);

  S=(double *)malloc(Rank*sizeof(double));
  assert(S != NULL);
/*
  Copy data
*/
  cblas_zcopy(N*Rank,aU,1,U,1);
  cblas_zcopy(M*Rank,aV,1,V,1);
/*
  LAPACK preparations
*/
  NL=(int )N;
  ML=(int )M;
  RankL=(int )Rank;
  LWork=(int )10*lmax(N,M);

//  UL=(MKL_Complex16 *)U;
//  VL=(MKL_Complex16 *)V;
//  TauUL=(MKL_Complex16 *)TauU;
//  TauVL=(MKL_Complex16 *)TauV;
//  WorkL=(MKL_Complex16 *)Work;
//  RL=(MKL_Complex16 *)R;
//  URL=(MKL_Complex16 *)UR;
//  VRL=(MKL_Complex16 *)VR;
//
//  aUL=(MKL_Complex16 *)aU;
//  aVL=(MKL_Complex16 *)aV;
  UL=(double complex *)U;
  VL=(double complex *)V;
  TauUL=(double complex *)TauU;
  TauVL=(double complex *)TauV;
  WorkL=(double complex *)Work;
  RL=(double complex *)R;
  URL=(double complex *)UR;
  VRL=(double complex *)VR;

  aUL=(double complex *)aU;
  aVL=(double complex *)aV;
/*
  QR-decomposition of the matrix U
*/
  zgeqrf_(&NL,&RankL,UL,&NL,TauUL,WorkL,&LWork,&Info);
  assert(Info == 0);
/*
  QR-decomposition of the matrix V
*/
  zgeqrf_(&ML,&RankL,VL,&ML,TauVL,WorkL,&LWork,&Info);
  assert(Info == 0);
/*
  Form the matrix R
*/
  zset(Rank*Rank,0.0+0.0*I,R,1);

  for (j=0; j < Rank; j++)
    for (i=0; i < Rank; i++)
      for (k=lmax(i,j); k < Rank; k++)
        R[j*Rank+i]+=U[k*N+i]*conj(V[k*M+j]);
/*
  SVD-decomposition of the matrix R
*/
  zgesvd_("S","S",&RankL,&RankL,RL,&RankL,S,URL,&RankL,VRL,&RankL,WorkL,&LWork,RWork,&Info);
  assert(Info == 0);
/*
  Frobenius norm of the matrix R and criterion
*/
  FroNorm=0.0;

  for (i=0; i < Rank; i++)
    FroNorm+=S[i]*S[i];

  Crit=0.25*Eps*Eps*FroNorm;
/*
  Compression check
*/
  Sum=0.0;
  i=Rank-1;

  do
  {
    NewRank=i+1;
    Sum+=S[i]*S[i];
    i-=1;
  }
  while ((Sum <= Crit) && (NewRank >= 0));
/*
  Successful compression ?
*/
  if (NewRank < Rank)
  {
/*
  Form orthogonal matrix U
*/
    for (j=0; j < Rank; j++)
    {
      cblas_zcopy(Rank,UR+j*Rank,1,aU+j*N,1);
      zset(N-Rank,0.0+0.0*I,aU+j*N+Rank,1);
    }
    zunmqr_("L","N",&NL,&RankL,&RankL,UL,&NL,TauUL,aUL,&NL,WorkL,&LWork,&Info);
    assert(Info == 0);
/*
  Form orthogonal matrix V
*/
    for (j=0; j < Rank; j++)
      for (i=0; i < Rank; i++)
        VR[j*Rank+i]=conj(VR[j*Rank+i]);

    for (j=0; j < Rank; j++)
    {
      cblas_zcopy(Rank,VR+j*Rank,1,aV+j,M);
      zset(M-Rank,0.0+0.0*I,aV+j*M+Rank,1);
    }

    zunmqr_("L","N",&ML,&RankL,&RankL,VL,&ML,TauVL,aVL,&ML,WorkL,&LWork,&Info);
    assert(Info == 0);
/*
  Scale the matrix U with the singular values
*/
    for (j=0; j < Rank; j++)
      for (i=0; i < N; i++)
        aU[j*N+i]*=S[j]+0.0*I;
/*
  New rank
*/
    *aRank=NewRank;
  }
/*
  Free memory
*/
  free(U); free(V); free(TauU); free(TauV); free(Work); 
  free(RWork); free(R); free(UR); free(VR); free(S);
}
