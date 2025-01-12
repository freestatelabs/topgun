#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define ERRMAX 1e-12
#define ITMAX 1000
const double pi =  M_PI;

double ellipK(double k2) {

    if (fabs(k2 - 1.0) <= ERRMAX) return 1e99;
    if (fabs(k2 - 0.0) <= ERRMAX) return 0.5*pi;

    int it = 1; 
    double err = 2*ERRMAX; 
    double a0 = 1.0; 
    double g0 = sqrt(1 - k2);
    double a1, g1;

    while ((it < ITMAX) && (err > ERRMAX)) {
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0);
        a0 = a1; g0 = g1;
        err = fabs(a0 - g0); it += 1;
    }
    // printf("\tExiting with it = %i\n", it);
    return pi/(2*a0);
}

double ellipKfast(double k2) {

    if (fabs(k2 - 1.0) <= ERRMAX) return 1e99;
    if (fabs(k2 - 0.0) <= ERRMAX) return 0.5*pi;

    int it = 1; 
    double err = 2*ERRMAX; 
    double a0 = 1.0; 
    double g0 = sqrt(1 - k2);
    double a1, g1;

    // 3 iterations 
    a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
    a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
    a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;

    err = fabs(a0 - g0);

    while ((err > ERRMAX) && it < ITMAX) {
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0);
        a0 = a1; g0 = g1;
        err = fabs(a0 - g0); it += 1;
    }

    return pi/(2*a0);
}

const int Nmax = 1e9;
const int N = 1e7;         // number of data points on the range
const double k2min = 0;
const double k2max = 1.0-1e-10;
const int Nit = Nmax/N;

int main() {
    clock_t start; 

    double* k2 = aligned_alloc(sizeof(double), sizeof(double)*N);
    double* K = aligned_alloc(sizeof(double), sizeof(double)*N);
    // Set k2 to be (0,1)
    double range = k2max - k2min;
    for (int i=0; i<N; i++) {
        k2[i] = i * range/(N - 1);
    }

    start = clock();
    for (int j=0; j<Nit; j++)
        for (int i=0; i<N; i++) {
            K[i] = ellipK(k2[i]);
            //printf("K[%.3f] = %.3f\n", k2[i], K[i]);
        }
    double ref_time = (double)1000*(clock() - start)/CLOCKS_PER_SEC/Nit;   // ms
    printf("Ref elapsed time =  %.3f ms\n", ref_time);

    start = clock();
    for (int j=0; j<Nit; j++)
        for (int i=0; i<N; i++) {
            K[i] = ellipKfast(k2[i]);
            // printf("K[%.3f] = %.3f\n", k2[i], K[i]);
        }

    double fast_time = (double)1000*(clock() - start)/CLOCKS_PER_SEC/Nit;
    printf("Fast elapsed time = %.3f ms\n", fast_time);
    printf("Speedup = %.2f %", 100*(ref_time - fast_time)/ref_time);
    return 1;
}