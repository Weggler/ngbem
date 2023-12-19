/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                19.09.2010                                   *
\*****************************************************************************/
/*
  Last change 19.09.2010
*/
/*
  Includes for all functions 
*/

# include "Includes.h"

void WriteSparseMatrix(char *FileName,SparseCOOMatrix SMat)
{
/*
  Local variables
*/
  long INZ;
  FILE *ResultFile;

  ResultFile=fopen(FileName,"w");

  fprintf(ResultFile,"%8ld %8ld %8ld\n",SMat.NRow,SMat.NColumn,SMat.NumNZ);
 
  for (INZ=0; INZ < SMat.NumNZ; INZ++)
    fprintf(ResultFile,"%8ld %8ld %25.15e %25.15e\n",SMat.Rows[INZ],SMat.Columns[INZ],creal(SMat.Values[INZ]),cimag(SMat.Values[INZ]));

  fclose(ResultFile);
}

 
