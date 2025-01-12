#include "time.h"
#include <stdio.h>
#include <stdlib.h>

const int Ni = 1000;
const int N = 48;


void matmul(const float *a, const float *_b, float *c, int n) {
    float *b = calloc(n*n, sizeof(float));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b[i * n + j] = _b[j * n + i];
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                c[i * n + j] += a[i * n + k] * b[j * n + k]; // <- note the indices
            }
        }
    }
    free(b);
}


int main() {
    float total_time = 0.0;

    float* a = calloc(N*N, sizeof(float));
    float* b = calloc(N*N, sizeof(float));
    float* c = calloc(N*N, sizeof(float));

    for (int i=0; i<N; i++){
        a[i] = (float)rand();
        b[i] = (float)rand();
        c[i] = (float)rand();
    }

    for (int i = 0; i < Ni; i++)
    {
        clock_t start = clock();
        matmul(a, b, c, 48);
        total_time += (float)(clock() - start)/CLOCKS_PER_SEC;
    }

    printf("Avg time per matmul: %.3f ms\n", total_time/Ni*1e3);
    
    free(a);
    free(b);
    free(c);
    return 0;
}