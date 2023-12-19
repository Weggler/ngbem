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

long ReadGeometry(Node *Nodes, Triangle *Elements)
{
    long i,Dummy,GN,N1,N2,N3;
    char String[1024];
    FILE *GeometryFile;
    
    GeometryFile=fopen(SurfaceFile,"r");
/*
  Search for the number of nodes
*/
    do
	fgets(String,sizeof(String),GeometryFile);
    while(strstr(String,"BEG_NODL_DATA") == 0);
/*
  Number of nodes
*/
    fgets(String,sizeof(String),GeometryFile);
/*
  Read nodes
*/
    for (i=0; i < NumNodes; i++)
    {
	fscanf(GeometryFile,"%ld %lf %lf %lf",&GN,
               &Nodes[i].X[0],&Nodes[i].X[1],&Nodes[i].X[2]);
	Nodes[i].GlNum=GN-1;
    }
/*
  Search for number of elements
*/
    do
	fgets(String,sizeof(String),GeometryFile);
    while(strstr(String,"BEG_ELEM_DATA") == 0);
/*
  Number of elements
*/
    fgets(String,sizeof(String),GeometryFile);

    for (i=0; i < NumElements; i++)
    {   
        fgets(String,sizeof(String),GeometryFile);
	sscanf(String,"%ld %ld %ld %ld %ld %ld",&GN,
               &Dummy,&Dummy,&N1,&N2,&N3);
	Elements[i].GlNum=GN-1;
        Elements[i].Nodes[0]=Nodes+N1-1;
	Elements[i].Nodes[1]=Nodes+N2-1;
	Elements[i].Nodes[2]=Nodes+N3-1;
    } 

    fclose(GeometryFile);

    return 0;
}
