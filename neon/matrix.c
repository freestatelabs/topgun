#include <stdlib.h> 
#include <stdio.h>
#include "matrix.h"

Matrix *mzeros(size_t m, size_t n) {
    Matrix *matrix = malloc(sizeof(Matrix)); 
    double *values = calloc(m*n, sizeof(double));
    matrix->m = m; 
    matrix->n = n; 
    matrix->values = values; 
    return matrix; 
} 

Matrix *mfroma(const double values[], size_t m, size_t n) {
    Matrix *matrix = mzeros(m, n); 
    for (size_t i=0; i<m; i++) {
        for (size_t j=0; j<n; j++) {
            matrix->values[i*n + j] = values[i*n+j];
        }
    }
    return matrix;
}


void mprint(Matrix *matrix) {
    printf("Matrix of size %li x %li:\n", matrix->m, matrix->n);
    for (size_t i=0; i<matrix->m; i++) {
        printf("\t");
        for (size_t j=0; j<matrix->n; j++) {
            printf("%f ", matrix->values[i*matrix->n+j]);
        }
        printf("\n");
    }
}

int mfree(Matrix *matrix) {
    free(matrix->values);
    free(matrix);
    return 0;
};