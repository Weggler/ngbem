/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                13.10.2011                                   *
\*****************************************************************************/
/*
  Last change 21.10.2011
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void JoinSets(long NA,long *A,long *NB,long **B)
{
  /*
    Local variables
  */
  long i,j;

//  printf("JoinSets: %ld %ld\n", NA, *NB);
//  fflush(stdout);

  assert((NA > 0) && (*NB >= 0));

  if(*NB == 0)
    *B=(long *)malloc(NA*sizeof(long));
  else
    *B=(long *)realloc(*B,((*NB)+NA)*sizeof(long));
  assert(*B != NULL);

  for (i=0; i < NA; i++)
    {
      j=0;

      while ((j < *NB)) 
	{
	  if ((*B)[j] == A[i]) break; 
          j++;
	}

      if (j == *NB)
	{
	  (*B)[*NB]=A[i];
	  (*NB)++;
	}
    }

  *B=(long *)realloc(*B,(*NB)*sizeof(long));
  assert(*B != NULL);
}
