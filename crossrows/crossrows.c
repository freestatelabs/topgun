#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int allocate_matrix_calloc(float** a, float** b, float** c, int N);
int allocate_matrix_malloc(float** a, float** b, float** c, int N);

int main(int argc, char* argv[]) {

    if (argc > 1) {
        for (int i=0; i<argc; i++) {
            printf("Arg #%i: %s\n", i, argv[i]);
        }
    }
    else {
        printf("No command line arguments.\n");
    }

    int N = 1000000;
    int Nit = 10000; 
    struct timespec start, end;
    float* a, *b, *c; 
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i=0; i<Nit; i++) {

        int error_code = allocate_matrix_malloc(&a, &b, &c, N); 
        a[0] = 1.0; 
        b[0] = 2.0; 
        c[0] = 3.0;

    }
    free(a);
    free(b); 
    free(c);
    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) +
                     (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("Calloc time: %.3e\n", elapsed);

    printf("Success!\n");
    return 0;
}

int allocate_matrix_calloc(float** a, float** b, float** c, int N) {
    *a = calloc(N, sizeof(float));
    *b = calloc(N, sizeof(float));
    *c = calloc(N, sizeof(float));

    if ((!*a) || (!*b) || (!*c)) {
        return -1;
    }
    else {
        return 0;
    }
}

int allocate_matrix_malloc(float** a, float** b, float** c, int N) {
    *a = malloc(N*sizeof(float));
    *b = malloc(N*sizeof(float));
    *c = malloc(N*sizeof(float));

    if ((!*a) || (!*b) || (!*c)) {
        return -1;
    }
    else {
        return 0;
    }
}