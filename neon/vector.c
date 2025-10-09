#include <stdlib.h>
#include <stdio.h>
#include "vector.h"

Vector *vzeros(size_t n) {
    Vector *vector = malloc(sizeof(Vector)); 
    double *values = calloc(n, sizeof(double));
    vector->n = n; 
    vector->values = values; 
    return vector;
} 

Vector *vfroma(double values[], size_t n) {
    Vector *vector = vzeros(n); 
    for (size_t i=0; i<n; i++) {
        vector->values[i] = values[i];
    }
    return vector;
}

void vprint(Vector *vector) {
    printf("Vector of length %li:\n", vector->n);

    for (size_t i=0; i<vector->n; i++) {
        printf("\t%f\n", vector->values[i]);
    }
}

int vfree(Vector *vector) {
    free(vector->values);
    free(vector);
    return 0;
}