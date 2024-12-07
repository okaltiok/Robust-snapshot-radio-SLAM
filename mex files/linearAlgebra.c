/*
 * linearAlgebra.c
*/

#include "linearAlgebra.h"
#include "math.h"

double pi_to_pi(double angle) {
    angle = fmod(angle, 2*M_PI);
    if (angle > M_PI)
    {
        angle -= 2*M_PI;
    }
    else if (angle < -M_PI)
    {
        angle += 2*M_PI;
    }
    return angle;
}

double rotation_matrix(double theta, double* R) {
    double c, s;
    
    c = cos(theta);
    s = sin(theta);
    
    R[0] = c; 
    R[1] = s;
    R[2] = -s;
    R[3] = c; 
}


bool choleskyFactorization(double *L, double *A, int dim)
{   
    double tmp;
    for (int i = 0; i < dim*dim; i++) 
        L[i] = 0.0;
    
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < j; k++)
                sum += L[i+k*dim] * L[j+k*dim];

            if (i == j)
            {
                tmp = A[i+i*dim] - sum;
                if (tmp <= 0)
                    return false;
                
                L[i+j*dim] = sqrt(tmp);
            }
            else
            {
                L[i+j*dim] = (A[i+j*dim] - sum) / L[j+j*dim] ;   
            }
        }
    }
    return true;
}

void matrixMultiply(int rA, int cA, int rB, int cB, int transpose, int opt, double *A, double *B, double *C)
{
    if (transpose == 0)
    {
        for (int i = 0; i < rA; i++) {
            for (int j = 0; j < cB; j++) {
                double sum = 0.0;
                for (int k = 0; k < cA; k++)
                {
                    sum = sum + A[i + rA * k] * B[j * rB + k];
                }
                switch(opt){
                    case 0:
                        C[i + rA * j] = sum;
                        break;
                    case 1:
                        C[i + rA * j] += sum;
                        break;
                    case -1:
                        C[i + rA * j] -= sum;
                        break;
                }
            }
        }
    }  
    else if (transpose == 1)
    {
        for (int i = 0; i < rA; i++) {
            for (int j = 0; j < rB; j++) {
                double sum = 0.0;
                for (int k = 0; k < cA; k++)
                {
                    sum = sum + A[i + rA * k] * B[j + rB * k];   
                }
                switch(opt){
                    case 0:
                        C[i + rA * j] = sum;
                        break;
                    case 1:
                        C[i + rA * j] += sum;
                        break;
                    case -1:
                        C[i + rA * j] -= sum;
                        break;
                }
            }
        }
    }
    else if (transpose == 2)
    {
        for (int i = 0; i < cA; i++) {
            for (int j = 0; j < cB; j++) {
                double sum = 0.0;
                for (int k = 0; k < rA; k++)
                {
                    sum = sum + A[i * rA + k] * B[j * rB + k];
                }
                switch(opt){
                    case 0:
                        C[i + cA * j] = sum;
                        break;
                    case 1:
                        C[i + cA * j] += sum;
                        break;
                    case -1:
                        C[i + cA * j] -= sum;
                        break;
                }
            }
        }
    }  
}

bool matrixInverse(double *invA, double *detA, int *iPivot, double *A, int N)
{   
    int i, j, k, imax;
    double maxA, ptr, absA;
    
    // Initialize arrays
    for (i = 0; i < N; i++)
        iPivot[i] = i;
    iPivot[N] = 0;
    
    for (int i = 0; i < N*N; i++)
        invA[i] = 0.0;

    
    // Compute LU decomposition
    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;
        
        for (k = i; k < N; k++){
            if ((absA = fabs(A[k + i*N])) > maxA) { 
                maxA = absA;
                imax = k;
            }
        }
        
        if (imax != i) {
            //pivoting iPivot
            j = iPivot[i];
            iPivot[i] = iPivot[imax];
            iPivot[imax] = j;

            //pivoting rows of A
            for (k = 0; k < N; k++){
                ptr = A[i + k*N];
                A[i + k*N] = A[imax + k*N];
                A[imax + k*N] = ptr;
            }
            //counting pivots starting from 0 (for determinant)
            iPivot[N]++;
        }
        
        for (j = i + 1; j < N; j++) {
            A[j + i*N] /= A[i + i*N];
            
            for (k = i + 1; k < N; k++)
                A[j + k*N] -= A[j + i*N] * A[i + k*N];
        }
    }
    
    
    // Compute determinant
    *detA = A[0];
    for (int i = 1; i < N; i++)
        *detA *= A[i + i*N];
    
    if (*detA == 0)
        return false;
    *detA =  iPivot[N] % 2 == 0 ? *detA : -*detA;
    
    
    // Compute matrix inverse
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            invA[i + j*N] = iPivot[i] == j ? 1.0 : 0.0;
            
            for (int k = 0; k < i; k++)
                invA[i + j*N] -= A[i + k*N] * invA[k + j*N];
        }
        
        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                invA[i + j*N] -= A[i + k*N] * invA[k + j*N];
            
            invA[i + j*N] /= A[i + i*N];
        }
    }
    
    return true;
}
