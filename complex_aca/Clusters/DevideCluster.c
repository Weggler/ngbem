//------------------------------------------------------------------
//                                                                  
//  routine name  -  DevideCluster
//                                                                  
//------------------------------------------------------------------
//                                                                  
//  last revision -  Aug 12
//
//  purpose       -             
//                                                                  
//  in:                                                             
//             G  - cluster weights
//             X  - cluster points
//     NClusters  - number of the new clusters
//          IClu  - number of actual cluster
//      Clusters  - cluster tree = list of Clusters
//        Permu   - Permutation of the points in cluster wrt global
//                  ordering in X,G (?)
//       ParamACA - parameters driving the devision process  
//                                                                  
//  out: 
//                                                                  
//------------------------------------------------------------------

/* Includes for all functions */
# include "../CInclude/Includes.h"

void DevideCluster(double *G,double *X,long *aNClusters,long IClu,Cluster **aClusters, long *Permu,ACAParameters ParamACA)
{
/* 
  Local variables 
*/
  int iprint = 0;
  long i,j,NClu,ClBegin,ClEnd,Num1,Num2,ICheck,IEnd,NClusters,IClu1,IClu2;
  long jMin;
  double Mom[10],XMin[3],XMax[3],Centre[3],Diag[3],Diff[3],XT[3];
  double DiagLength,Crit,alpha,Radius,Dist,DistMin;
  Cluster *Clusters;
/* 
  LAPACK variables 
*/
  int N,LDA,LWORK,INFO;
  double A[9],W[3],WORK[8];
/* 
  Number of elements 
*/
  Clusters =*aClusters;
  NClusters=*aNClusters;
/* 
  Number of elements 
*/
  NClu   =Clusters[IClu].Number;
  ClBegin=Clusters[IClu].PermuPos;
  ClEnd  =Clusters[IClu].PermuPos+NClu-1;
  if(iprint == 1)
    printf("DevideCluster NClu=%ld, ClBegin=%ld, ClEnd=%ld\n", NClu,ClBegin,ClEnd);
/* 
  Moments 
*/
  dset(10,0.0,Mom,1);

  for(i=ClBegin; i <= ClEnd; i++)
  {
    j=Permu[i];

    Mom[0]+=G[j];
    Mom[1]+=G[j]*X[3*j];
    Mom[2]+=G[j]*X[3*j+1];
    Mom[3]+=G[j]*X[3*j+2];
    Mom[4]+=G[j]*X[3*j]*X[3*j];
    Mom[5]+=G[j]*X[3*j]*X[3*j+1];
    Mom[6]+=G[j]*X[3*j]*X[3*j+2];
    Mom[7]+=G[j]*X[3*j+1]*X[3*j+1];
    Mom[8]+=G[j]*X[3*j+1]*X[3*j+2];
    Mom[9]+=G[j]*X[3*j+2]*X[3*j+2];
  }
/* 
  Centre 
*/
  cblas_dcopy(3,Mom+1,1,Centre,1);
  alpha=1.0/Mom[0];
  cblas_dscal(3,alpha,Centre,1);
/* 
  Covariance matrix 
*/
  A[0]=Mom[4]-Centre[0]*Mom[1];
  A[1]=A[3]=Mom[5]-Centre[0]*Mom[2];
  A[2]=A[6]=Mom[6]-Centre[0]*Mom[3];
  A[4]=Mom[7]-Centre[1]*Mom[2];
  A[5]=A[7]=Mom[8]-Centre[1]*Mom[3];
  A[8]=Mom[9]-Centre[2]*Mom[3];
/* 
  LAPACK preparations 
*/
  N=3;
  LDA=3;
  LWORK=8;
/* 
  LAPACK for eigevalues and eigenvectors 
*/
  dsyev_("V","U",&N,A,&LDA,W,WORK,&LWORK,&INFO);
  assert(INFO == 0);
/* 
  Radius and bounding box 
*/
  Radius=0.0;
   
  for(i=ClBegin; i <= ClEnd; i++)
  {
    j=Permu[i];
/* 
  X-Centre 
*/
    cblas_dcopy(3,X+3*j,1,Diff,1);
    cblas_daxpy(3,-1.0,Centre,1,Diff,1);
/* 
  Search for the basis element in the cluster 
*/
    Dist=cblas_dnrm2(3,Diff,1);
   
    if(i == ClBegin)
    {
      jMin=j;
      DistMin=Dist;
    }
    else if (Dist < DistMin)
    {
      jMin=j;
      DistMin=Dist;
    }
/* 
  X in local coordinates 
*/
    XT[0]=cblas_ddot(3,A,1,Diff,1);
    XT[1]=cblas_ddot(3,A+3,1,Diff,1);
    XT[2]=cblas_ddot(3,A+6,1,Diff,1);
    
    Radius=max(Radius,cblas_dnrm2(3,XT,1));
    
    if (i == ClBegin)
	{
	  XMin[0]=XT[0];
	  XMax[0]=XT[0];

	  XMin[1]=XT[1];
	  XMax[1]=XT[1];

	  XMin[2]=XT[2];
	  XMax[2]=XT[2];
	}
    else 
	{
	  XMin[0]=min(XMin[0],XT[0]);
	  XMax[0]=max(XMax[0],XT[0]);

	  XMin[1]=min(XMin[1],XT[1]);
	  XMax[1]=max(XMax[1],XT[1]);

	  XMin[2]=min(XMin[2],XT[2]);
	  XMax[2]=max(XMax[2],XT[2]);
	}
  }
/* 
  Diagonal 
*/
  cblas_dcopy(3,XMax,1,Diag,1);
  cblas_daxpy(3,-1.0,XMin,1,Diag,1);
  DiagLength=cblas_dnrm2(3,Diag,1);
  if (DiagLength > 0.0)
  {
    alpha=1.0/DiagLength;
    cblas_dscal(3,alpha,Diag,1);
  }
/* 
  Copy geometrical information 
*/
  Clusters[IClu].BasisElement=jMin;

  Clusters[IClu].Radius=Radius;
  cblas_dcopy(3,W,1,Clusters[IClu].EVal,1);
  cblas_dcopy(9,A,1,Clusters[IClu].EVec,1);
  cblas_dcopy(3,XMin,1,Clusters[IClu].XMin,1);
  cblas_dcopy(3,XMax,1,Clusters[IClu].XMax,1);
  cblas_dcopy(3,Centre,1,Clusters[IClu].Centre,1);
  Clusters[IClu].DiagLength=DiagLength;
/* 
  Devision 
*/
  if (NClu <= ParamACA.NCluMin)
  {
/* 
  Small cluster, no devision 
*/
    Clusters[IClu].Son1=-1;
    Clusters[IClu].Son2=-1;
    return;
  }
  else 
  {
/* 
  Big cluster, devision 
*/
    Crit=cblas_ddot(3,Centre,1,A+6,1);
   
    Num1=0;
    Num2=0;
   
    ICheck=ClBegin;
    IEnd=ClEnd;
   
    for (i=0; i < NClu; i++)
    {
      j=Permu[ICheck];
    
      if (cblas_ddot(3,X+3*j,1,A+6,1) >= Crit)
      {
        Num1+=1;
        ICheck+=1;
      }
      else
      {
        Num2+=1;
        Permu[ICheck]=Permu[IEnd];
        Permu[IEnd]=j;
        IEnd-=1;
      }
    }
/* 
  Non-separable cluster 
*/
    if ((Num1 == 0) || (Num2 == 0))
	{
	  Clusters[IClu].Son1=-1;
	  Clusters[IClu].Son2=-1;
	  return;
	} 
    else
    {
/* 
  New clusters 
*/
      NClusters+=2;
              
      Clusters=realloc(Clusters,NClusters*sizeof(Cluster));
      assert(Clusters != NULL);
/* 
  Information and further devision 
*/
      IClu1=NClusters-2;
      IClu2=NClusters-1;
/* 
  Save the first son 
*/
      Clusters[IClu].Son1=IClu1;
      
      Clusters[IClu1].Level=Clusters[IClu].Level+1;            
      Clusters[IClu1].Father=IClu;            
      Clusters[IClu1].Number=Num1;            
      Clusters[IClu1].PermuPos=ClBegin;            
/* 
  Save the second son 
*/
      Clusters[IClu].Son2=IClu2;
           
      Clusters[IClu2].Level=Clusters[IClu].Level+1;            
      Clusters[IClu2].Father=IClu;            
      Clusters[IClu2].Number=Num2;            
      Clusters[IClu2].PermuPos=ClBegin+Num1;            
    } 
  }

  *aClusters =Clusters;
  *aNClusters=NClusters;
}
