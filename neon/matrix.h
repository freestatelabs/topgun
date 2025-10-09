/*  Private header file for Matrix operations
*/

#ifndef MATRIX_H 

#include <stdlib.h>

typedef struct {
    size_t m; 
    size_t n; 
    double *values;
} Matrix;

#endif 