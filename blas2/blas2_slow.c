
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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
    double **values;
} Matrix; 

Vector *vrand(size_t nrows, double min, double max) {
    Vector *vector = malloc(sizeof(Vector));
    double *values = calloc(nrows, sizeof(double));
    for (size_t i=0; i<nrows; i++) {
        values[i] = drand(min, max);
    }
    vector->nrows = nrows; 
    vector->values = values;
    return vector;
}

Matrix *mrand(size_t nrows, size_t ncols, double min, double max) {
    Matrix *matrix = malloc(sizeof(Matrix));
    double **values = malloc(nrows*sizeof(double*));
    for (size_t i=0; i<nrows; i++) {
        values[i] = calloc(ncols, sizeof(double));
        for (size_t j=0; j<ncols; j++) {
            values[i][j] = drand(min, max);
        }
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
    for (size_t i=0; i<matrix->nrows; i++) {
        free(matrix->values[i]);
    }
    free(matrix->values);
    free(matrix);
}

// A*x = b
int mmultv(Matrix *A, Vector *x, Vector *b) {
    if (A->ncols != x->nrows) { return -1; }
    for (size_t i=0; i<A->nrows; i++) {
        double acc = 0.0;
        for (size_t j=0; j<A->ncols; j++) {
            acc += A->values[i][j]*x->values[j];
            //printf("i, j = %li, %li\n", i, j);
        }
        b->values[i] = acc;
    }
    return 0;
}

int main() {
    printf("Testing blas2 slow...\n");

    size_t n = 10000;
    size_t niterations = 100;
    Matrix *A = mrand(n, n, -100, 100);
    Vector *x = vrand(n, -100, 100);
    Vector *b = vrand(n, 0, 0);

    clock_t start = clock();
    int test = 0;
    double checksum = 0;
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

