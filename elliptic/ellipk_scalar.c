#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define ERRMAX 1e-12
#define ITMAX 1000
const double pi =  M_PI;

// Return the average of an array of integers
float average(int* Nt, int N) {
    float avg = 0;
    for (int i=0; i<N; i++) avg += Nt[i];
    return avg/N;
}

// Return the minimum of an array of integers
int min_i(int* Nt, int N) {
    int m = ITMAX*2; 
    for (int i=0; i<N; i++) {
        if (Nt[i] < m) m = Nt[i]; 
    }
    return m;
}

double min_d(double* Nt, int N) {
    double m = ITMAX*2; 
    for (int i=0; i<N; i++) {
        if (Nt[i] < m) m = Nt[i]; 
    }
    return m;
}

// Return the maximum of an array of integers
int max_i(int* Nt, int N) {
    int m = 0; 
    for (int i=0; i<N; i++) {
        if (Nt[i] > m) m = Nt[i]; 
    }
    return m;
}

double ellipK_ref(double k2, int* Nt) {

    if (fabs(k2 - 1.0) <= ERRMAX) return 1e99;
    if (fabs(k2 - 0.0) <= ERRMAX) return 0.5*pi;

    int it = 0; 
    double err = 2*ERRMAX; 
    double a0 = 1.0; 
    double g0 = sqrt(1 - k2);
    double a1, g1;

    while ((it < ITMAX) && (err > ERRMAX)) {
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0);
        a0 = a1; g0 = g1;
        err = fabs(a0 - g0); it += 1;
    }
    
    *Nt = it;
    return pi/(2*a0);
}

double ellipKfast(const double k2) {

    //double k2m1 = k2 - 1.0;
    // if ((k2m1 <= ERRMAX) {
    //     if(k2m1 > 0)) return 1e99;
    // if (k2 < 0.0 || k2 > 1.0) return -1;
    // if (k2 <= ERRMAX) return 0.5*pi;

    // double err = 2*ERRMAX; 
    double a0 = 1.0; 
    double g0 = sqrt(1 - k2);
    double a1, g1;

    // Precompute iterations using magic numbers derived from testing
    if (k2 < 0.0094) {
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
    }
    else if (k2 < 0.326) {
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
    }
    else if (k2 < 0.931)
    {
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
    }
    else if (k2 < 0.999) { 
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
    }
    else {        
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
        a1 = 0.5*(a0+g0); g1 = sqrt(a0*g0); a0 = a1; g0 = g1;
    }

    //double err = fabs(a0 - g0);
    //if (err > ERRMAX) printf("ERROR at k2 = %f; err = %e\n", k2, err);
    return pi/(2*a0);
}


const int N = 1e2;         // number of data points on the range
const double k2min = 0.0;
const double k2max = 1.0-ERRMAX;
const int Nit = 1e0;

int main() {
    clock_t start; 

    double* k2 = aligned_alloc(sizeof(double), sizeof(double)*N);
    double* K = aligned_alloc(sizeof(double), sizeof(double)*N);
    double* Kfast = aligned_alloc(sizeof(double), sizeof(double)*N);
    int* Nt = calloc(sizeof(int), sizeof(int)*N); 
    // Set k2 to be (0,1)
    double range = k2max - k2min;
    for (int i=0; i<N; i++) {
        k2[i] = k2min + i * range/(N - 1);
    }

    start = clock();
    for (int j=0; j<Nit; j++)
        for (int i=0; i<N; i++) {
            K[i] = ellipK_ref(k2[i], &Nt[i]);
            //printf("K[%.12f] = %.3f; iterations = %i\n", k2[i], K[i], Nt[i]);
        }
    double ref_time = (double)1000*(clock() - start)/CLOCKS_PER_SEC/Nit;   // ms
    printf("Ref elapsed time =  %.3f ms\n", ref_time);
    printf("Average iterations: %.1f\n", average(Nt, N));
    printf("Max iterations:     %i\n", max(Nt, N));
    printf("Min iterations:     %i\n", min(Nt, N)); 

    start = clock();
    for (int j=0; j<Nit; j++)
        for (int i=0; i<N; i++) {
            Kfast[i] = ellipKfast(k2[i]);
            // printf("K[%.12f] = %.3f, Kfast[%.12f] = %.3f; iterations = %i\n", k2[i], K[i], k2[i], Kfast[i], Nt[i]);
        }

    // Compute difference 
    double* diff = aligned_alloc(sizeof(double), N*sizeof(double));
    for (int i=0; i<N; i++) {
        diff[i] = K[i] - Kfast[i];
    }
    double avgdiff = average(diff, N);
    double maxdiff = max(diff, N);
    free(diff);


    double fast_time = (double)1000*(clock() - start)/CLOCKS_PER_SEC/Nit;
    printf("Fast elapsed time = %.3f ms\n", fast_time);
    printf("Speedup = %.2f %", 100*(ref_time - fast_time)/ref_time);
    printf("Max/avg diff: %.12f/%.12f\n", avgdiff, maxdiff);

    free(k2); free(K); free(Kfast); free(Nt);
    return 1;
}