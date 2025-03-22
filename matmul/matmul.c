#include "time.h"
#include <stdio.h>
#include <stdlib.h>

const int Ni = 1000;
const int m = 100; 
const int n = 100; 
const int p = 100; 


void matmul(const float *a, const float *b, float *c, int m, int n, int p) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            for (int k = 0; k < n; k++)
                c[j * n + i] += a[k * n + i] * b[j * n + k];
}

void matmul2(const float *_a, const float *b, float *c, int m, int n, int p) {
    // Transposed version
    float* a = calloc(m*n, sizeof(float)); 

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            a[j*p + i] = _a[i*p + j];
        }
    }

    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            for (int k = 0; k < n; k++)
                c[j * n + i] += a[k * n + i] * b[k * n + j];

    free(a);
}

int main() {
    // float total_time = 0.0;

    float* a = calloc(m*n, sizeof(float));
    float* b = calloc(n*p, sizeof(float));
    float* c = calloc(m*p, sizeof(float));

    for (int i=0; i<m*n; i++){
        a[i] = (float)rand()/RAND_MAX;
    }
    for (int i=0; i<n*p; i++){
        b[i] = (float)rand()/RAND_MAX;
    }


    struct timespec start, end;
    float total_time = 0.0;


    for (int i = 0; i < Ni; i++)
    {
        for (int i=0; i<m*p; i++) { 
            c[i] = 0.0;
        }
        // clock_t start = clock();
        // matmul(a, b, c, m, n, p);
        // total_time += (float)(clock() - start)/CLOCKS_PER_SEC;
        clock_gettime(CLOCK_MONOTONIC, &start);
        matmul2(a, b, c, m, n, p);
        clock_gettime(CLOCK_MONOTONIC, &end);
        total_time  += (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) * 1e-9;
    }

    printf("Avg time per matmul: %.3f us\n", (total_time/Ni)*1e6);
    
    free(a);
    free(b);
    free(c);
    return 0;
}