/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                03.04.2009                                   *
\*****************************************************************************/
/*
  Includes for the main function 
*/

# include "IncludesMain.h"

void ReadPointSet(long *N,double **X,double **G)
{
/*
  Local variables
*/
  long i,Dummy;
  char String[4096];
  FILE *GeometryFile;

  GeometryFile=fopen(SurfaceFile,"r");
/*
  Number of points
*/
  fgets(String,sizeof(String),GeometryFile);
  sscanf(String,"%ld",N);

  (*X)=(double *)malloc(3*(*N)*sizeof(double));
  assert((*X) != NULL);

  (*G)=(double *)malloc((*N)*sizeof(double));
  assert((*G) != NULL);
/*
  Read points and weights
*/
  for (i=0; i < (*N); i++)
    fscanf(GeometryFile,"%ld %lf %lf %lf %lf",
	   &Dummy,(*X)+3*i,(*X)+3*i+1,(*X)+3*i+2,(*G)+i);
}
