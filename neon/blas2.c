/*  BLAS L2 routines; in this case, just various implementations of DGEMV:
        (D)ouble-precision (GE)neral (M)atrix-(V)ector multiplication

    Variants: 
    * Naive (complete)
    * Blocked (in-progress)
    * Manual SIMD, not blocked (todo)
    * Manual SIMD, blocked (todo)
*/
#include <stdlib.h>

#include "matrix.h"
#include "vector.h"

int dgemv_naive(const Matrix *restrict A, const Vector *restrict x, Vector *restrict y) {

    // x must have same number of columns as A; y must have same number of rows
    if ((A->n != x->n) || (A->m != y->n)) { return -1; }

    // Outer loop across rows of A 
    for (size_t i=0; i<A->m; i++) {
        double acc = 0.0;
        // Inner loop over elements 
        for (size_t j=0; j<A->n; j++) {
            acc += A->values[i*A->n+j]*x->values[j];
        }
        y->values[i] = acc;
    }
    return 0;
}