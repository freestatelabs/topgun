/* GEMV (matrix-vector multiplication) using manual neon intrinsics
*/

#include <stdlib.h> 
#include <stdio.h>
#include <arm_neon.h>

#include "linalg.h"

int main() {
    const size_t m = 4; 
    const size_t n = 3; 
    double values[n] = {1, 2, 3};
    Vector *x = vfroma(values, n);
    Vector *y = vzeros(m);

    double mvalues[m*n] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; 
    Matrix *A = mfroma(mvalues, m, n);
    int success = dgemv_naive(A, x, y); 
    vprint(y);

    mfree(A);
    vfree(x);
    vfree(y);
    return 0;
}
