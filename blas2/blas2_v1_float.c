/*  v1 
    Faster matrix-vector multiply. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Accelerate/Accelerate.h>

float drand(float min, float max) {
    return min + (max - min) * rand() / RAND_MAX;
}

typedef struct {
    size_t nrows;
    float *values; 
} Vector;

typedef struct {
    size_t nrows; 
    size_t ncols; 
    float *values;
} Matrix; 

Vector *vrand(size_t nrows, float min, float max) {
    Vector *vector = malloc(sizeof(Vector));
    float *values = calloc(nrows, sizeof(float));
    for (size_t i=0; i<nrows; i++) {
        values[i] = drand(min, max);
    }
    vector->nrows = nrows; 
    vector->values = values;
    return vector;
}

Matrix *mrand(size_t nrows, size_t ncols, float min, float max) {
    Matrix *matrix = malloc(sizeof(Matrix));
    float *values = malloc(nrows*ncols*sizeof(float));
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
        float acc = 0.0;
        for (size_t j=0; j<A->ncols; j++) {
            acc += A->values[i*A->ncols + j]*x->values[j];
        }
        b->values[i] = acc;
    }
    return 0;
}

int mmultv_accelerate(const Matrix *A,
                       const Vector *x,
                             Vector *b)
{
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

int main() {
    printf("Testing blas2 v1 (float)...\n");
    // cblas_set_num_threads(1);

    size_t n = 10000;
    size_t niterations = 100;
    Matrix *A = mrand(n, n, -100, 100);
    Vector *x = vrand(n, -100, 100);
    Vector *b = vrand(n, 0, 0);

    clock_t start = clock();
    int test = 0;
    float checksum = 0;
    for (size_t i=0; i<niterations; i++) {
        test += mmultv(A, x, b); 
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