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
int min(int* Nt, int N) {
    int m = ITMAX*2; 
    for (int i=0; i<N; i++) {
        if (Nt[i] < m) m = Nt[i]; 
    }
    return m;
}

// Return the maximum of an array of integers
int max(int* Nt, int N) {
    int m = 0; 
    for (int i=0; i<N; i++) {
        if (Nt[i] > m) m = Nt[i]; 
    }
    return m;
}

/*
    double ellipK_ref(double k2)

Reference implementation on scalar inputs. Saves the number of iterations to Nt
*/
double ellipK_ref(double k2, int* Nt) {

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

    *Nt = it;
    return pi/(2*a0);
}


const int N = 1e6;         // number of data points on the range
const double k2min = 0;
const double k2max = 1.0 - ERRMAX;
const int Nit = 1000;

int main() {
    clock_t start; 

    double* k2 = aligned_alloc(sizeof(double), sizeof(double)*N);
    double* K = aligned_alloc(sizeof(double), sizeof(double)*N);
    int* Nt = calloc(sizeof(int), sizeof(int)*N); 

    // Set k2 to be (0,1)
    double range = k2max - k2min;
    for (int i=0; i<N; i++) {
        k2[i] = i * range/(N - 1);
    }

    start = clock();
    for (int j=0; j<Nit; j++)
        for (int i=0; i<N; i++) {
            K[i] = ellipK_ref(k2[i], &Nt[i]);
            //printf("K[%.3f] = %.3f; iterations = %i\n", k2[i], K[i], Nt[i]);
        }
    double ref_time = (double)1000*(clock() - start)/CLOCKS_PER_SEC/Nit;   // ms
    printf("Ref elapsed time =  %.3f ms\n", ref_time);
    printf("Average iterations: %.1f\n", average(Nt, N));
    printf("Max iterations:     %i\n", max(Nt, N));
    printf("Min iterations:     %i\n", min(Nt, N)); 

    free(k2); free(K); free(Nt);
    return 1;
}
