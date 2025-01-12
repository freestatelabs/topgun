#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int Ni = 1000;
const int N = 48;
const int B = 8;

typedef float vec __attribute__ ((vector_size(32))); // 8 floats per vector

vec* alloc(int n) {
    vec* ptr = (vec*) aligned_alloc(32, sizeof(vec) * n); // Aligned allocation for vectors
    memset(ptr, 0, sizeof(vec) * n); // Initialize to 0
    return ptr;
}

float hsum(vec s) {
    float res = 0;
    for (int i = 0; i < B; i++) {
        res += s[i];
    }
    return res;
}

void matmul(const float* _a, const float* _b, float* c, int n) {
    int nB = (n + B - 1) / B; // Number of 8-element vectors per row (rounded up)

    vec* a = alloc(n * nB);
    vec* b = alloc(n * nB);

    // Move both matrices to the aligned region
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i * nB + j / B][j % B] = _a[i * n + j]; // Store in a 32-byte aligned vector
            b[j * nB + i / B][i % B] = _b[j * n + i]; // Store in b with transposition
        }
    }

    // Matrix multiplication
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            vec s = {0};

            // Vertical summation
            for (int k = 0; k < nB; k++) {
                s += a[i * nB + k] * b[j * nB + k];
            }

            c[i * n + j] += hsum(s);
        }
    }

    free(a);
    free(b);
}

int main() {
    srand(time(NULL)); // Seed the random number generator

    float total_time = 0.0;
    float* a = calloc(N * N, sizeof(float)); // Allocate memory for matrices a, b, and c
    float* b = calloc(N * N, sizeof(float));
    float* c = calloc(N * N, sizeof(float));

    for (int i = 0; i < N * N; i++) {
        a[i] = (float)rand();
        b[i] = (float)rand();
        c[i] = (float)rand();
    }

    for (int i = 0; i < Ni; i++) {
        clock_t start = clock();
        matmul(a, b, c, N); // Pass N as the matrix dimension
        total_time += (float)(clock() - start) / CLOCKS_PER_SEC;
    }

    printf("Avg time per matmul: %.3f ms\n", total_time / Ni * 1e3);

    free(a);
    free(b);
    free(c);

    return 0;
}
