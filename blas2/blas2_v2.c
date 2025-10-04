/*  v2
    Faster matrix-vector multiply. Use memalign
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
// #include <omp.h>

double drand(double min, double max) {
    return min + (max - min) * rand() / RAND_MAX;
}

typedef struct {
    size_t nrows;
    double *values; 
} Vector;

typedef struct {
    size_t nrows; 
    size_t ncols; 
    double *values;
} Matrix; 

#define ALIGNMENT 1024
Vector *vrand(size_t nrows, double min, double max) {
    Vector *vector = malloc(sizeof(Vector));
    double *values = calloc(nrows, sizeof(double));
    // int success = posix_memalign((void**)&values, ALIGNMENT, nrows*sizeof(double));
    // printf("Success? %i\n", success);
    for (size_t i=0; i<nrows; i++) {
        values[i] = drand(min, max);
    }
    vector->nrows = nrows; 
    vector->values = values;
    return vector;
}

Matrix *mrand(size_t nrows, size_t ncols, double min, double max) {
    Matrix *matrix = malloc(sizeof(Matrix));
    double *values = malloc(nrows*ncols*sizeof(double));    
    //int success = posix_memalign((void**)&values, ALIGNMENT, ncols*nrows*sizeof(double));
    // printf("Success? %i\n", success);
    for (size_t i=0; i<nrows*ncols; i++) {
        values[i]= drand(min, max);
    }

    matrix->nrows = nrows; 
    matrix->ncols = ncols; 
    matrix->values = values; 
    return matrix;
}

void vfree(Vector *vector) {
    free(vector->values);
    free(vector);
}

void mfree(Matrix *matrix) {
    free(matrix->values);
    free(matrix);
}

// A*x = b
int mmultv(Matrix *restrict A, Vector *restrict x, Vector *restrict b) {
    if (A->ncols != x->nrows) { return -1; }
    for (size_t i=0; i<A->nrows; i++) {
        double acc = 0.0;
        for (size_t j=0; j<A->ncols; j++) {
            acc += A->values[i*A->ncols + j]*x->values[j];
        }
        b->values[i] = acc;
    }
    return 0;
}

bool vequal(Vector *v1, Vector *v2) {

    if (v1->nrows != v2->nrows) { return false; }

    double err = 0.0; 


    for (size_t i=0; i<v1->nrows; i++) {
        err += fabs(v1->values[i] - v2->values[i]);
    }

    printf("Err = %e\n", err);

    if (err > 1e-6) {
        return false;
    }
    else {
        return true;
    }
}

int mmultv_accelerate(const Matrix *A,
                       const Vector *x,
                             Vector *b)
{
    
    // BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);
    /* Accelerate expects column‑major storage.  Either:
       • Allocate A in column‑major order, or
       • Pass the transposed flag (CblasNoTrans) and swap the leading dimension.
       Here we keep the row‑major layout and tell BLAS to treat it as transposed. */
    cblas_dgemv(CblasRowMajor,          // layout matches our storage
                CblasNoTrans,           // we really want y = A*x
                (int)A->nrows,
                (int)A->ncols,
                1.0,                     // alpha
                A->values,
                (int)A->ncols,          // lda = stride between rows
                x->values,
                1,                       // incX
                0.0,                     // beta
                b->values,
                1);                      // incY
    return 0;
}


int mmultv_tile(const Matrix *A, const Vector *x, Vector *y, 
                size_t tile_m, size_t tile_n) {

    if ((A->ncols != x->nrows) || (x->nrows != y->nrows)) { return -1; }

    // Outer loop down each row of the matrix
    for (size_t i0 = 0; i0 < A->nrows; i0 += tile_m) {

        // Determine the number of rows to work over
        // Make sure we only work to the end of the matrix
        size_t imax = (i0 + tile_m > A->nrows) ? A->nrows : i0 + tile_m;

        for (size_t j0 = 0; j0 < A->ncols; j0 += tile_n) {
            size_t jmax = (j0 + tile_n > A->ncols) ? A->ncols : j0 + tile_n;

            /* Inner kernel works on a TILE_M × TILE_N block */
            for (size_t i = i0; i < imax; ++i) {
                double acc = y->values[i];           // start with existing y[i]

                const double *row = &A->values[i * A->ncols + j0];
                for (size_t j = j0; j < jmax; ++j) {
                    acc += row[j - j0] * x->values[j];
                }
                y->values[i] = acc;
            }
        }
    }
    return 0;
}

int main() {
    printf("Testing blas2 v2...\n");
    // cblas_set_num_threads(1);

    size_t n = 10000;          // standard test has been 100k
    size_t niterations = 100;
    size_t tile_m = 8;
    size_t tile_n = 8;

    Matrix *A = mrand(n, n, -100, 100);
    Vector *x = vrand(n, -100, 100);
    Vector *b = vrand(n, 0, 0);
    Vector *b2 = vrand(n, 0, 0);
    printf("all allocs succeeded\n");

    // Test 
    mmultv(A, x, b); 
    mmultv_tile(A, x, b2, tile_m, tile_n); 
    // mmultv_accelerate(A, x, b2);
    bool is_equal = vequal(b, b2); 
    if (is_equal) {
        printf("Test passed!\n");
    }
    else {
        printf("Test failed!\n");
    }


    clock_t start = clock();
    int test = 0;
    double checksum = 0;
    for (size_t i=0; i<niterations; i++) {
        //test += mmultv(A, x, b); 
        test += mmultv_tile(A, x, b, tile_m, tile_n);
        // test += mmultv_accelerate(A, x, b);
        checksum += b->values[2];
    }
    printf("checksum = %f\n", checksum);
    clock_t end = clock(); 
    double elapsed = (double)(end - start)/CLOCKS_PER_SEC;
    elapsed /= (double)niterations; 

    printf("Time per multiply: %.3f ms\n", 1e3*elapsed);
    printf("Test: %i\n", test);     // make sure this is 0
    mfree(A);
    vfree(x);
    vfree(b);
}

/*
16.4ms for cblas/accelerate

33 ms for julia (BLAS), 1 thread
19 ms for julia (BLAS), 2 threads
15 ms for julia (BLAS), 4 threads
13 ms for julia (BLAS), 6 threads

123ms for naive implementation

tile timings: (N, M) = TIME
(8, 8) = 30 ms
(16, 8) = 34 ms
(4, 4) = 42 ms
(2, 2) = 83 ms
(4, 8) = 42 ms
(8, 12) = 29.5 ms
(10, 10) = 35.5 ms
(12, 12) = 31.1 ms
*/


