/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                          Laplace equation in 3D                             *
*                                                                             *
*                                12.11.2006                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "Includes.h"

long ElementInfo(Triangle *Triangles)
{
    double U12[3],U23[3],U31[3],Q,d12,d23,d31,d,h,ttau,tstar;
    long i;
    double *X1,*X2,*X3;
    Triangle *TR;

    for(i=0; i < NumElements; i++)
    {
/*
  Actual tirangle
*/
	TR=Triangles+i;
/*
  Adresses of the edges of the triangle
*/
	X1=TR->Nodes[0]->X;
	X2=TR->Nodes[1]->X;
	X3=TR->Nodes[2]->X;
/*
  Sides of the triangle
*/
	cblas_dcopy(3,X2,1,U12,1);
	cblas_daxpy(3,-1.0,X1,1,U12,1);
        d12=cblas_dnrm2(3,U12,1);

	cblas_dcopy(3,X3,1,U23,1);
	cblas_daxpy(3,-1.0,X2,1,U23,1);
        d23=cblas_dnrm2(3,U23,1);

	cblas_dcopy(3,X1,1,U31,1);
	cblas_daxpy(3,-1.0,X3,1,U31,1);
        d31=cblas_dnrm2(3,U31,1);

	d=max(max(d12,d23),d31);
 /*
  Unit outer normal vector
*/
        dvecpr(U31,U12,TR->Norm);
        Q=cblas_dnrm2(3,TR->Norm,1);
        if (d*d >= Q*X_TLOSS)
        {
            printf("  Degenerated triangle %d \n",(int)(i+1));
            fflush(stdout);
	    return i+1;
        }
        else
        {
/*
  Scaling normal vector
*/
            cblas_dscal(3,1.0/Q,TR->Norm,1);
/*
  Area of the tirangle
*/
	    TR->Area=0.5*Q;
/*
  Local mesh size
*/
            h=sqrt(Q);
/*
  Max-Min mesh size, uniformity
*/
            if (i == 0) 
            {
		HMax=h;
		HMin=h;
		CB=d/h;
            }
	    else
            {
		HMax=max(HMax,h);
		HMin=min(HMin,h);
		CB=max(CB,d/h);
	    }
/*
  Local system of coordinates in X1
*/
            ttau=cblas_dnrm2(3,U23,1);

	    cblas_dcopy(3,U23,1,TR->LocCS[0][1],1);
	    cblas_dscal(3,1.0/ttau,TR->LocCS[0][1],1);

            tstar=-cblas_ddot(3,U12,1,TR->LocCS[0][1],1);

            cblas_dcopy(3,X2,1,TR->LocCS[0][0],1);
	    cblas_daxpy(3,tstar,TR->LocCS[0][1],1,TR->LocCS[0][0],1);
            cblas_daxpy(3,-1.0,X1,1,TR->LocCS[0][0],1);

            TR->Height[0]=cblas_dnrm2(3,TR->LocCS[0][0],1);
            
	    cblas_dscal(3,1.0/TR->Height[0],TR->LocCS[0][0],1);

            TR->TanAngle[0][0]=-tstar/TR->Height[0];
            TR->TanAngle[0][1]=(ttau-tstar)/TR->Height[0];
/*
  Local system of coordinates in X2
*/
            ttau=cblas_dnrm2(3,U31,1);

	    cblas_dcopy(3,U31,1,TR->LocCS[1][1],1);
	    cblas_dscal(3,1.0/ttau,TR->LocCS[1][1],1);

            tstar=-cblas_ddot(3,U23,1,TR->LocCS[1][1],1);

            cblas_dcopy(3,X3,1,TR->LocCS[1][0],1);
	    cblas_daxpy(3,tstar,TR->LocCS[1][1],1,TR->LocCS[1][0],1);
            cblas_daxpy(3,-1.0,X2,1,TR->LocCS[1][0],1);

            TR->Height[1]=cblas_dnrm2(3,TR->LocCS[1][0],1);
            
	    cblas_dscal(3,1.0/TR->Height[1],TR->LocCS[1][0],1);

            TR->TanAngle[1][0]=-tstar/TR->Height[1];
            TR->TanAngle[1][1]=(ttau-tstar)/TR->Height[1];
/*
  Local system of coordinates in X3
*/
            ttau=cblas_dnrm2(3,U12,1);

	    cblas_dcopy(3,U12,1,TR->LocCS[2][1],1);
	    cblas_dscal(3,1.0/ttau,TR->LocCS[2][1],1);

            tstar=-cblas_ddot(3,U31,1,TR->LocCS[2][1],1);

            cblas_dcopy(3,X1,1,TR->LocCS[2][0],1);
	    cblas_daxpy(3,tstar,TR->LocCS[2][1],1,TR->LocCS[2][0],1);
            cblas_daxpy(3,-1.0,X3,1,TR->LocCS[2][0],1);

            TR->Height[2]=cblas_dnrm2(3,TR->LocCS[2][0],1);
            
	    cblas_dscal(3,1.0/TR->Height[2],TR->LocCS[2][0],1);

            TR->TanAngle[2][0]=-tstar/TR->Height[2];
            TR->TanAngle[2][1]=(ttau-tstar)/TR->Height[2];
	}
    }
    return 0;
}
