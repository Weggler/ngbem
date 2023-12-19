//                                                                  
// S. Rjasanow: Adaptive Cross Approximation                 
//
//------------------------------------------------------------------
//                                                                  
//  header name  -  Cluster
//                                                                  
//------------------------------------------------------------------
//                                                                  
//  last revision -  Mar 12
//
//   purpose      - definition of struct Cluster and struct Pari
//
//------------------------------------------------------------------

// Level          - cluster level
// Father,Son1,2  - to reconstruct the tree
// Number         - dimension (?)
// PermuPos       - to reconstruct global position
//                                                                  
// Radius         - geometrical
// EVal[3]        - data
// EVec[9]        - for 
// XMin[3]        - its 
// XMax[3]        - construction (?)
// Centre[3]      -
// DiagLength     -
//
typedef struct sCluster
{
    long Level;
    long Father;
    long Son1;
    long Son2;
    long Number;
    long PermuPos;

    double Radius;
    double EVal[3];
    double EVec[9];
    double XMin[3];
    double XMax[3];
    double Centre[3];
    double DiagLength;
} Cluster;

//  Clu1,2 - the cluster numbers of the pair
//  Type   - its type, i.e.,
typedef enum {Admissible=0,Dense=1,Hierarchical=2} BlockType;
typedef struct sPair
{
    long Clu1;
    long Clu2;
    BlockType Type;
} Pair;
