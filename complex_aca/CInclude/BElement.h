/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
**                                                                             *
*                                20.12.2008                                   *
\*****************************************************************************/
/*
  Definitions of Nodes, Edges and Triangles
*/

struct sCluster;

typedef struct sNode
{
    long GlNum;
    double X[3];
    long NumNb;
    struct sTriangle **Nb;
    long  *KNb;
    void *BC;
} Node;

typedef struct sEdge
{
    long GlNum;
    Node *Nodes[2];
    double Length;
    long NumNb;
    struct sTriangle **Nb;
    void *BC;
} Edge;

typedef struct sTriangle
{
    long GlNum;
    Node *Nodes[3];
    Edge *Edges[3];
    double Norm[3];
    double Area;
    double Height[3];
    double LocCS[3][2][3];
    double TanAngle[3][2];
    struct sTriangle *Nb[3];
    void *BC;
} Triangle;
