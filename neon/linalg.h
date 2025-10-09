/*  Public interface for `neon` linear algebra library
    (c) 2025 ryan@freestatelabss
*/

#ifndef LINALG_H 

#include <stdlib.h>

// a 2D row-major Matrix of doubles
typedef struct {} Matrix; 
// a 1D vector of doubles
typedef struct {} Vector; 

Vector *vzeros(size_t n); 
Vector *vfroma(double values[], size_t n);
void vprint(Vector *vector);
int vfree(Vector *vector);

Matrix *mzeros(size_t m, size_t n); 
Matrix *mfroma(double values[], size_t m, size_t n);
void mprint(Matrix *matrix);
int mfree(Matrix *matrix);

// Matrix-vector multiplication; A*x = y
int dgemv_naive(const Matrix *restrict A, const Vector *restrict x, Vector *restrict y);

#endif