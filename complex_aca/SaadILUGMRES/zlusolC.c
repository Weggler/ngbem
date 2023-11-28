#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "zheads.h"
#include "zprotos.h"
#include "zdefs.h"

int zlusolC( complex double *y, complex double *x, iluptr lu )
{
/*----------------------------------------------------------------------
 *    performs a forward followed by a backward solve
 *    for LU matrix as produced by iluk
 *    y  = right-hand-side
 *    x  = solution on return
 *    lu = LU matrix as produced by iluk.
 *--------------------------------------------------------------------*/
    int n = lu->n, i, j, nzcount, *ja;
    complex double *D;
    csptr L, U;

    L = lu->L;
    U = lu->U;
    D = lu->D;

   /* Block L solve */
    for( i = 0; i < n; i++ ) {
        x[i] = y[i];
        nzcount = L->nzcount[i];
        ja = L->ja[i];
        for( j = 0; j < nzcount; j++ ) {
            x[i] -= x[ja[j]] * L->ma[i][j];
        }
    }

    /* Block -- U solve */
    for( i = n-1; i >= 0; i-- ) {
        nzcount = U->nzcount[i];
        ja = U->ja[i];
        for( j = 0; j < nzcount; j++ ) {
            x[i] -= x[ja[j]] * U->ma[i][j];
        }
        x[i] *= D[i];
    }
    return (0);
}
