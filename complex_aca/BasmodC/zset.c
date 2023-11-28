# include "../CInclude/Includes.h"

void zset(long N,double complex A,double complex *X, long IX)
{ 
    long i, j;

    if(IX==1)
	for(i=0; i<N; i++)
	    X[i]=A;
    else
	for(j=IX*N, i=0; i<j; i+=IX)
	    X[i]=A;

    return;
}
