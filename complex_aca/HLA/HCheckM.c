/*****************************************************************************\
*                                                                             *
*                          S. Rjasanow C-Software                             *
*                                                                             *
*                       Adaptive Cross Approximation                          *
*                                                                             *
*                                07.01.2011                                   *
\*****************************************************************************/
/*
  Last change 25.11.2010
*/
/*
  Includes for all functions 
*/

# include "../CInclude/Includes.h"

void HCheckM(CBlock *A)
{
  double *Err;

  if((*A).Type  == Null)
    {
      if(((*A).A11 != NULL) ||
	  ((*A).A21 != NULL) ||
 	  ((*A).A12 != NULL) ||
	  ((*A).A22 != NULL) ||
	  ((*A).Data != NULL) ||
	  ((*A).Rank != 0) ||
	  ((*A).UData != NULL) ||
	  ((*A).VData != NULL))
	{
	  printf("In HCheckM : Wrong N-Matrix\n");
	  fflush(stdout);
	  printf("      Matrix type : %ld\n",(long)(*A).Type);
	  fflush(stdout);
	  printf("      A11,A21,A12,A22 = %p %p %p %p\n",(*A).A11,(*A).A21,(*A).A12,(*A).A22);
	  fflush(stdout);
	  printf("      Data = %p\n",(*A).Data);
	  fflush(stdout);
	  printf("      Rank, UData,VData = %ld %p %p\n",(*A).Rank,(*A).UData,(*A).VData);
	  fflush(stdout);
/*
  Stop
*/
	  Err=NULL;
          (*Err)=1.0;
	}
    }
  else if((*A).Type  == Hierarchical)
    {
      if(
         ((*A).A11 == NULL) ||
	 ((*A).A21 == NULL) ||
 	 ((*A).A12 == NULL) ||
	 ((*A).A22 == NULL) ||
	 ((*A).Data != NULL) ||
	 ((*A).Rank != 0) ||
	 ((*A).UData != NULL) ||
	 ((*A).VData != NULL)
        )
	{
	  printf("In HCheckM : Wrong H-Matrix\n");
	  fflush(stdout);
	  printf("      Matrix type : %ld\n",(long)(*A).Type);
	  fflush(stdout);
	  printf("      A11,A21,A12,A22 = %p %p %p %p\n",(*A).A11,(*A).A21,(*A).A12,(*A).A22);
	  fflush(stdout);
	  printf("      Data = %p\n",(*A).Data);
	  fflush(stdout);
	  printf("      Rank, UData,VData = %ld %p %p\n",(*A).Rank,(*A).UData,(*A).VData);
	  fflush(stdout);
/*
  Stop
*/
	  Err=NULL;
          (*Err)=1.0;
	}
    }
  else if ((*A).Type == Dense)
    {
      if(((*A).A11 != NULL) ||
	  ((*A).A21 != NULL) ||
 	  ((*A).A12 != NULL) ||
	  ((*A).A22 != NULL) ||
	  ((*A).Data == NULL) ||
	  ((*A).Rank != 0) ||
	  ((*A).UData != NULL) ||
	  ((*A).VData != NULL))
	{
	  printf("In HCheckM : Wrong D-Matrix\n");
	  fflush(stdout);
	  printf("      Matrix type : %ld\n",(long)(*A).Type);
	  fflush(stdout);
	  printf("      A11,A21,A12,A22 = %p %p %p %p\n",(*A).A11,(*A).A21,(*A).A12,(*A).A22);
	  fflush(stdout);
	  printf("      Data = %p\n",(*A).Data);
	  fflush(stdout);
	  printf("      Rank, UData,VData = %ld %p %p\n",(*A).Rank,(*A).UData,(*A).VData);
	  fflush(stdout);
/*
  Stop
*/
	  Err=NULL;
          (*Err)=1.0;
	}
    }
  else if ((*A).Type == Admissible)
    {
      if(((*A).A11 != NULL) ||
	  ((*A).A21 != NULL) ||
 	  ((*A).A12 != NULL) ||
	  ((*A).A22 != NULL) ||
	  ((*A).Data != NULL) ||
	  ((*A).Rank == 0) ||
	  ((*A).UData == NULL) ||
	  ((*A).VData == NULL))
	{
	  printf("In HCheckM : Wrong A-Matrix\n");
	  fflush(stdout);
	  printf("      Matrix type : %ld\n",(long)(*A).Type);
	  fflush(stdout);
	  printf("      A11,A21,A12,A22 = %p %p %p %p\n",(*A).A11,(*A).A21,(*A).A12,(*A).A22);
	  fflush(stdout);
	  printf("      Data = %p\n",(*A).Data);
	  fflush(stdout);
	  printf("      Rank, UData,VData = %ld %p %p\n",(*A).Rank,(*A).UData,(*A).VData);
	  fflush(stdout);
/*
  Stop
*/
	  Err=NULL;
          (*Err)=1.0;
	}
    }
}
