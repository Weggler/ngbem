/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                          Laplace equation in 3D                             *
*                                                                             *
*                                15.12.2003                                   *
\*****************************************************************************/
/*
  Includes for all functions 
*/

# include "Includes.h"

void FindNodeNeighbours(Node *Nodes, Triangle *Triangles)
{
  long i,j,NumNbAll;
/*
  Count the neighbours 
*/
  for(i=0; i < NumNodes; i++)
    Nodes[i].NumNb=0;

  for(i=0; i < NumElements; i++)
    {
      Nodes[Triangles[i].Nodes[0]->GlNum].NumNb++;
      Nodes[Triangles[i].Nodes[1]->GlNum].NumNb++;
      Nodes[Triangles[i].Nodes[2]->GlNum].NumNb++;
    }
/*
  Number of neighbours
*/
  NumNbAll=0;
  for(i=0; i < NumNodes; i++)
    NumNbAll+=Nodes[i].NumNb;
/*
  Allocation of the memory for the neighbours for each node
*/
  Nodes[0].Nb=malloc(NumNbAll*sizeof(Triangle*));
  assert(Nodes[0].Nb != NULL);
/*
  Allocation of the memory for the neighbours for each node
*/
  Nodes[0].KNb=malloc(NumNbAll*sizeof(long));
  assert(Nodes[0].KNb != NULL);
/*
  Fill the field of neighbours
*/
  for(i=1; i < NumNodes; i++)
    {
      Nodes[i].Nb=Nodes[i-1].Nb+Nodes[i-1].NumNb;
      Nodes[i].KNb=Nodes[i-1].KNb+Nodes[i-1].NumNb;
      Nodes[i-1].NumNb=0;
    }
  Nodes[NumNodes-1].NumNb=0;
 
  for(i=0; i < NumElements; i++)
    { 
      j=Triangles[i].Nodes[0]->GlNum;
      Nodes[j].Nb[Nodes[j].NumNb]=Triangles+i;
      Nodes[j].KNb[Nodes[j].NumNb++]=0;

      j=Triangles[i].Nodes[1]->GlNum;
      Nodes[j].Nb[Nodes[j].NumNb]=Triangles+i;
      Nodes[j].KNb[Nodes[j].NumNb++]=1;

      j=Triangles[i].Nodes[2]->GlNum;
      Nodes[j].Nb[Nodes[j].NumNb]=Triangles+i;
      Nodes[j].KNb[Nodes[j].NumNb++]=2;
    }
}
/*****************************************************************************\
*                                                                             *
*                                    The end                                  *
*                                                                             *
\*****************************************************************************/
