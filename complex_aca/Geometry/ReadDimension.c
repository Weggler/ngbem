/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                20.12.2008                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "Includes.h"

long ReadDimension(void)
{
    char String[4096];
    FILE *GeometryFile;
    
    GeometryFile=fopen(SurfaceFile,"r");
/*
  Search for the number of nodes
*/
    do
	fgets(String,sizeof(String),GeometryFile);
    while(strstr(String,"BEG_NODL_DATA") == 0);
/*
  Read the number of nodes
*/
    fgets(String,sizeof(String),GeometryFile);
    sscanf(String,"%ld",&NumNodes);
/*
  Search for number of elements
*/
    do
	fgets(String,sizeof(String),GeometryFile);
    while(strstr(String,"BEG_ELEM_DATA") == 0);
/*
  Read number of elements
*/
    fgets(String,sizeof(String),GeometryFile);
    sscanf(String,"%ld",&NumElements);

    fclose(GeometryFile);

    return 0;
}
