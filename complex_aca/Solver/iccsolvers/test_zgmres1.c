#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "zgmres1.h"

void mv(double complex *in,double complex *out,void *params) {
	int i;
	double complex ddd[10];
	for(i=0;i<10;i++) ddd[i]=1.0+1.0*I;
	for(i=0;i<10;i++) {
		out[i]=ddd[i]*in[i];
	}
}
void prn_func(int i, double err) {
	fprintf(stdout,"err[%4d]= %le\n",i,err);
	fflush(stdout);
}

int main() {
	int n=10;
	int i,itmax=100;
	double complex rhs[10], sol[10],tst[10];
	for(i=0;i<10;i++) rhs[i]=1.0+1.0*I;
	for(i=0;i<10;i++) sol[i]=0.0+0.0*I;
	i=zfgmres1(n, mv, NULL,  rhs, sol, 1.0e-10, 100, &itmax, prn_func);
	printf("\n%d %d\n\n",i,itmax);
	printf("solution\n");
	for(i=0;i<10;i++) fprintf(stdout,"%d %le %le\n",i,creal(sol[i]),cimag(sol[i]));
	mv(sol,tst,NULL);
	printf("             rhs                 A*x     \n");
	for(i=0;i<10;i++) fprintf(stdout,"%d %le %le %le %le \n",i,creal(rhs[i]),cimag(rhs[i]),creal(tst[i]),cimag(tst[i]));
	return 0;
}
