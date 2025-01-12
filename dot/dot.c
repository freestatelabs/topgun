#include "time.h"
#include <stdio.h>
#include <stdlib.h>
#include <arm_neon.h>

const int Ni = 1000; 
const int N = 100000;

// Convenience function to create a 2D array of zeros
float** zeros(int Nrows, int Ncols) {
    float** array; 
    array = calloc(Nrows, sizeof(float*));
    for (int i=0; i<Nrows; i++) {
        array[i] = calloc(Ncols, sizeof(float));
    }
    return array;
}

// De-allocate a 2D array
void free2darray(float** array, int Ncols) {
    for (int i=0; i<Ncols; i++) {
        free(array[i]);
    }
    free(array);
}


void dot(float** a, float** b, float* c){
    for (int i=0; i<N; i++)
    c[i] = a[i][0]*b[i][0] + a[i][1]*b[i][1] + a[i][2]*b[i][2];
}

void dot2(float* a1, float* a2, float* a3, float* b1, float* b2, float* b3, float* c) {
    for (int i=0; i<N; i++) {
        c[i] += a1[i] * b1[i];
    }
    for (int i=0; i<N; i++){
        c[i] += a2[i] * a2[i];
    }
    for (int i=0; i<N; i++) {
        c[i] += a3[i] * b3[i];
    }
}

// Neon-optimized dot product function
void dot3(float** a, float** b, float* c) {
    float32x4_t sum = vdupq_n_f32(0.0f); // Initialize a Neon register for sum

    for (int i = 0; i < N; i += 4) {
        // Load 4 elements from both arrays a and b into Neon registers
        float32x4_t a_values = vld1q_f32(&a[i][0]);
        float32x4_t b_values = vld1q_f32(&b[i][0]);

        // Perform element-wise multiplication
        float32x4_t mul = vmulq_f32(a_values, b_values);

        // Accumulate the results
        sum = vaddq_f32(sum, mul);
    }

    // Store the results into a normal array to reduce the sum (sum[0] + sum[1] + sum[2] + sum[3])
    float result[4];
    vst1q_f32(result, sum);
    c[0] = result[0] + result[1] + result[2] + result[3];
}

int main() {
    float total_time = 0.0;

    float** a = zeros(N,3);
    float** b = zeros(N,3);
    float* c = (float*)aligned_alloc(32, N*32);
    //float* a1 = calloc(N, sizeof(float));
    float* a1 = (float*)aligned_alloc(32, N*32);
    float* a2 = (float*)aligned_alloc(32, N*32);
    float* a3 = (float*)aligned_alloc(32, N*32);
    float* b1 = (float*)aligned_alloc(32, N*32);
    float* b2 = (float*)aligned_alloc(32, N*32);
    float* b3 = (float*)aligned_alloc(32, N*32);

    for (int i=0; i<N; i++){
        a1[i] = 1.0;
        a2[i] = 2.0;
        a3[i] = 3.0;
        b1[i] = 1.0;
        b2[i] = 2.0;
        b3[i] = 3.0;
        c[i] = 0.0;
    }

    for (int i = 0; i < Ni; i++)
    {
        clock_t start = clock();
        dot(a, b, c);
        total_time += (float)(clock() - start)/CLOCKS_PER_SEC;
    }

    printf("Avg time per dot: %.3f ms\n", total_time/Ni*1e3);
    
    total_time = 0.0;

    for (int i = 0; i < Ni; i++)
    {
        clock_t start = clock();
        dot2(a1, a2, a3, b1, b2, b3, c);
        total_time += (float)(clock() - start)/CLOCKS_PER_SEC;
    }

    printf("Avg time per dot2: %.3f ms\n", total_time/Ni*1e3);
    
    free2darray(a,3);
    free2darray(b,3);
    free(c);
    free(a1);
    free(a2);
    free(a3);
    free(b1);
    free(b2);
    free(b3);
    return 0;
}