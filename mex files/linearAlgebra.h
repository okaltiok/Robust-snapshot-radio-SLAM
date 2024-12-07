#ifndef LINEARALGEBRA_H_
#define LINEARALGEBRA_H_

#include <stdbool.h>

double pi_to_pi(double angle);

double rotation_matrix(double theta, double* R);

bool choleskyFactorization(double *L, double *A, int dim);

void matrixMultiply(int rA, int cA, int rB, int cB, int transpose, int opt, double *A, double *B, double *C);

// bool matrixInverse(double *invX, double *det, double *X, int dim);

bool matrixInverse(double *invA, double *detA, int *iPivot, double *A, int N);

#endif