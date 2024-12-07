//     This function computes the cost matrix, state estimates and outliers of a geometry-based snapshot SLAM solution
// 
//     Input:
//        permutations        - (n_k x M) permutation matrix
//        y                   - (3 x m_k) measurement matrix
//        tx                  - (3 x 1) TX state
//        theta               - (N x 1) RX heading
//        los_candidate       - (1 x 1) LOS index
//        gain                - (1 x m_k) channel gain
//        epsilon             - (1 x 1) error threshold for inliers
//
//     Output:
//        L                 - (M x N) containing the cost of the estimate
//        x_hat             - (3 x M x N) containing the state estimates
//        outliers          - (m_k x M x N) containing the outliers
//     
//     Author   : Ossi Kaltiokallio
//                Tampere University, Department of Electronics and
//                Communications Engineering
//                Korkeakoulunkatu 1, 33720 Tampere
//                ossi.kaltiokallio@tuni.fi
//     Last Rev : 15/11/2023
//     Tested   : '9.8.0.1359463 (R2020a) Update 1'
//     
//     Copyright notice: You are free to modify, extend and distribute 
//        this code granted that the author of the original code is 
//        mentioned as the original author of the code.

#include "config.h"            
#include "linearAlgebra.h"
#include "mex.h"
#include "string.h"
#include "math.h"
#include "stdbool.h"

void compute_parameters(double *tx, double* y, double theta, int los_candidate, double* gain, int m_k, double* MU, double* ETA, double* ETA_bar, double* HH, double* AA, double* bb);

void compute_model(double *x_hat, double *AA, double *bb, double *idx, int m_k);

bool is_feasible(double *x_hat, double *y, double *MU, double *ETA, double *HH, double *idx, int m_k, double eta2_threshold);

void compute_cost(double *nu2, double *x_hat, double *HH, double *MU, double *ETA_bar, int m_k);

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* inputs */
    double* permutations = mxGetDoubles(prhs[0]);
    double* y = mxGetDoubles(prhs[1]);
    double* tx = mxGetDoubles(prhs[2]);
    double* theta = mxGetDoubles(prhs[3]);
    int los_candidate = mxGetScalar(prhs[4]);
    double* gain = mxGetDoubles(prhs[5]);
    double epsilon = mxGetScalar(prhs[6]);
    double eta2_threshold = mxGetScalar(prhs[7]);
    
    /* dimension declarations */
    int M = (int) mxGetN(prhs[0]);    // number of permutations
    int n_k = (int) mxGetM(prhs[0]);  // number of elements in one permutation
    int N = (int) mxGetM(prhs[3]);    // number of RX headings
    int m_k = (int) mxGetN(prhs[1]);    // number of measurements
                
    /* create the output matrix */
    const mwSize ndim = 3;
    const mwSize dims1[3] = {x_dim,M,N};
    const mwSize dims2[3] = {m_k,M,N};

    plhs[0] = mxCreateDoubleMatrix(M,N, mxREAL);
    plhs[1] = mxCreateNumericArray(ndim, dims1, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(ndim, dims2, mxDOUBLE_CLASS, mxREAL);
    
    double* L = mxGetDoubles(plhs[0]);
    double* x_hat = mxGetDoubles(plhs[1]);
    double* outliers = mxGetDoubles(plhs[2]);
    
    /* help variables */
    mxArray* nu2_ptr = mxCreateDoubleMatrix(m_k, 1, mxREAL);
    mxArray* inliers_ptr = mxCreateDoubleMatrix(m_k, 1, mxREAL);
    mxArray* MU_ptr = mxCreateDoubleMatrix(p_dim, m_k, mxREAL);
    mxArray* ETA_ptr = mxCreateDoubleMatrix(p_dim, m_k, mxREAL);
    mxArray* ETA_bar_ptr = mxCreateDoubleMatrix(p_dim, m_k, mxREAL);
    mxArray* HH_ptr = mxCreateDoubleMatrix(p_dim, x_dim*m_k, mxREAL);
    mxArray* AA_ptr = mxCreateDoubleMatrix(x_dim, x_dim*m_k, mxREAL);
    mxArray* bb_ptr = mxCreateDoubleMatrix(x_dim, m_k, mxREAL);

    double* nu2 = mxGetDoubles(nu2_ptr);
    double* inliers = mxGetDoubles(inliers_ptr);
    double* MU = mxGetDoubles(MU_ptr);
    double* ETA = mxGetDoubles(ETA_ptr);
    double* ETA_bar = mxGetDoubles(ETA_bar_ptr);
    double* HH = mxGetDoubles(HH_ptr);
    double* AA = mxGetDoubles(AA_ptr);
    double* bb = mxGetDoubles(bb_ptr);

    
    double nan = mxGetNaN();
    int N_inliers;
    
    /* loop through all RX headings */
    for (int n=0; n < N; n++)
    {
                    
        /* compute parameters */         
        compute_parameters(tx, y, *theta, los_candidate, gain, m_k, MU, ETA, ETA_bar, HH, AA,  bb);
        
        /* loop through all permutations */
        for (int m=0; m < M; m++)
        {
            /* determine inliers */
            for (int i=0; i < m_k; i++)
                inliers[i] = 0;
            
            for (int i=0; i < n_k; i++)
                inliers[(int)(permutations[i]-1)] = 1;
            
            /* determine consensus set */
            compute_model(x_hat, AA, bb, inliers, m_k);
            N_inliers = 0;
            if (is_feasible(x_hat, y, MU, ETA, HH, inliers, m_k, eta2_threshold))
            {
                compute_cost(nu2, x_hat, HH, MU, ETA_bar, m_k);
                
                for (int i=0; i < m_k; i++)
                {
                    if (nu2[i] > epsilon)
                    {
                        outliers[i] = 1.0;
                        inliers[i] = 0;
                    }
                    else
                    {
                        inliers[i] = 1;
                        N_inliers++;
                    }
                }
            }
            
            /* fit model */
            *L = nan;
            if (N_inliers >= n_k)
            {
                compute_model(x_hat, AA, bb, inliers, m_k);
                
                if (is_feasible(x_hat, y, MU, ETA, HH, inliers, m_k, eta2_threshold))
                {
                    compute_cost(nu2, x_hat, HH, MU, ETA_bar, m_k);
                    
                    *L = 0.0;
                    for (int i=0; i < m_k; i++)
                    {
                        if (outliers[i] == 1)
                            *L += epsilon*gain[i];
                        else
                            *L += nu2[i]*gain[i];
                    }
                }
            }
            
            /* increment arrays */
            permutations += n_k;
            L ++;
            x_hat += x_dim;
            outliers += m_k;   
        }
        /* increment arrays */
        theta++;
        permutations -= n_k*M;
    }
    
    /* free arrays */
    mxDestroyArray(nu2_ptr);
    mxDestroyArray(inliers_ptr);
    mxDestroyArray(MU_ptr);
    mxDestroyArray(ETA_ptr);
    mxDestroyArray(ETA_bar_ptr);
    mxDestroyArray(HH_ptr);
    mxDestroyArray(AA_ptr);
    mxDestroyArray(bb_ptr);
}

void compute_parameters(double *tx, double* y, double theta, int los_candidate, double* gain, int m_k, double* MU, double* ETA, double* ETA_bar, double* HH, double* AA, double* bb)
{
    double rot_tx[DIM_P2];
    double rot_rx[DIM_P2];
    double tmp[DIM_P];
    double ui[DIM_P];
    double vi[DIM_P];
    double mu[DIM_P];
    double eta[DIM_P];
    double eta_bar[DIM_P];
    double H[DIM_PxX];
    double I[DIM_P2];
    double tmp2[DIM_PxX];
    double eta_bar2[DIM_P2];
    double A[DIM_X2];
    double b[DIM_X];
    double tau, aod, aoa;
    
    rotation_matrix(tx[2], rot_tx);
    rotation_matrix(theta, rot_rx);

    H[0] = 1.0;
    H[1] = 0.0;
    H[2] = 0.0;
    H[3] = 1.0;
    
    I[0] = 1.0;
    I[1] = 0.0;
    I[2] = 0.0;
    I[3] = 1.0;
    
    /* loop through measurements */
    for (int i=0; i < m_k; i++)
    {
        tau = y[i*h_dim];
        aod = y[i*h_dim+1];
        aoa = y[i*h_dim+2];
        
        tmp[0] = cos(aod);
        tmp[1] = sin(aod);
        matrixMultiply(p_dim, p_dim, p_dim, 1, 0, 0, rot_tx, tmp, ui);
        
        tmp[0] = cos(aoa);
        tmp[1] = sin(aoa);
        matrixMultiply(p_dim, p_dim, p_dim, 1, 0, 0, rot_rx, tmp, vi);
        
        mu[0] = tx[0] - tau*vi[0];
        mu[1] = tx[1] - tau*vi[1];
        
        H[4] = -vi[0];
        H[5] = -vi[1];
        
        if (i == los_candidate - 1)
        {
            eta[0] = 0.0;
            eta[1] = 0.0;
            eta_bar[0] = 0.0;
            eta_bar[1] = 0.0;
            
            matrixMultiply(p_dim, x_dim, p_dim, x_dim, 2, 0, H, H, A);
            matrixMultiply(p_dim, x_dim, p_dim, 1, 2, 0, H, mu, b);
        }
        else
        {
            eta[0] = ui[0] + vi[0];
            eta[1] = ui[1] + vi[1];
            eta_bar[0] = eta[0]/sqrt(pow(eta[0],2)+pow(eta[1],2));
            eta_bar[1] = eta[1]/sqrt(pow(eta[0],2)+pow(eta[1],2));
            
            memcpy(eta_bar2, I, SOD_P2);
            matrixMultiply(p_dim, 1, p_dim, 1, 1, -1, eta_bar, eta_bar, eta_bar2);
            matrixMultiply(p_dim, x_dim, p_dim, p_dim, 2, 0, H, eta_bar2, tmp2);
            
            matrixMultiply(x_dim, p_dim, p_dim, x_dim, 0, 0, tmp2, H, A);
            matrixMultiply(x_dim, p_dim, p_dim, 1, 0, 0, tmp2, mu, b);
        }
        
        for (int j=0; j < x_dim2; j++)
            A[j] *= gain[i];
        
        for (int j=0; j < x_dim; j++)
            b[j] *= gain[i];
        
        
        memcpy(MU+p_dim*i, mu, SOD_P);
        memcpy(ETA+p_dim*i, eta, SOD_P);
        memcpy(ETA_bar+p_dim*i, eta_bar, SOD_P);
        memcpy(HH+p_dim*x_dim*i, H, SOD_PxX);
        memcpy(AA+x_dim2*i, A, SOD_X2);
        memcpy(bb+x_dim*i, b, SOD_X);
    }
}

void compute_model(double *x_hat, double *AA, double *bb, double *idx, int m_k)
{
    double b[DIM_X];            // x_dim x 1
    double A[DIM_X2];           // x_dim x x_dim
    double Ainv[DIM_X2];        // x_dim x x_dim
    double Adet[DIM_ONE];       // 1 x 1
    int iPivot[DIM_X_PIVOT];    // (x_dim + 1) x 1
    int j;

    for (int ii=0; ii < x_dim; ii++)
        b[ii] = 0.0;

    for (int ii=0; ii < x_dim2; ii++)
        A[ii] = 0.0;

    /* loop through measurements */
    for (int i=0; i < m_k; i++)
    {
        if ((int)idx[i] == 1)
        {
            for (int ii=0; ii < x_dim; ii++)
                b[ii] += bb[i*x_dim+ii];
            
            for (int ii=0; ii < x_dim2; ii++)
                A[ii] += AA[i*x_dim2+ii];
        }
    }

    // invert A and compute determinant
    matrixInverse(Ainv, Adet, iPivot, A, x_dim);

    /* compute x = inv(A)*b */
    matrixMultiply(x_dim, x_dim, x_dim, 1, 0, 0, Ainv, b, x_hat);
}

bool is_feasible(double *x_hat, double *y, double *MU, double *ETA, double *HH, double *idx, int m_k, double eta2_threshold)
{
    int los_idx = 0;
    double min_tau = y[0];

    /* check if propagation distance is negative */
    for (int i=0; i < m_k; i++)
    {
        if (y[i*h_dim] < min_tau)
        {
            min_tau = y[i*h_dim];
            los_idx = i;
        }

        if (y[i*h_dim] - x_hat[2] < 0)
            return false;
    }   
    
    double eta2[DIM_ONE];            // x_dim x 1
    double tmp[DIM_P];            // x_dim x 1
    double gamma[DIM_ONE];
    
    /* loop through measurements */
    for (int i=0; i < m_k; i++)
    {
        if ((int)idx[i] == 1)
        {
            matrixMultiply(p_dim, 1, p_dim, 1, 2, 0, (ETA+i*p_dim), (ETA+i*p_dim), eta2);

            if (i == los_idx & eta2[0] < eta2_threshold)
            {
                gamma[0] = 1.0;
            }
            else
            {
                matrixMultiply(p_dim, x_dim, x_dim, 1, 0, 0, (HH+i*p_dim*x_dim), x_hat, tmp);

                for (int ii=0; ii < p_dim; ii++)
                    tmp[ii] -= MU[i*p_dim+ii];

                matrixMultiply(p_dim, 1, p_dim, 1, 2, 0, (ETA+i*p_dim), tmp, gamma);
                gamma[0] /= ((y[i*h_dim] - x_hat[2]) * eta2[0]);
            }

            if (gamma[0] < 0.0 | gamma[0] > 1.0)
            {
                return false;
            }
        }
    }
    return true;
}

void compute_cost(double *nu2, double *x_hat, double *HH, double *MU, double *ETA_bar, int m_k)
{
    double tmp[DIM_P];            // x_dim x 1
    double nu[DIM_P];           // x_dim x x_dim
    double coeff[DIM_ONE];            // x_dim x 1
    
    /* loop through measurements */
    for (int i=0; i < m_k; i++)
    {
        matrixMultiply(p_dim, x_dim, x_dim, 1, 0, 0, (HH+i*p_dim*x_dim), x_hat, tmp);
        
        for (int ii=0; ii < p_dim; ii++)
            tmp[ii] -= MU[i*p_dim+ii];
        
        matrixMultiply(p_dim, 1, p_dim, 1, 2, 0, (ETA_bar+i*p_dim), tmp, coeff);
        
        for (int ii=0; ii < p_dim; ii++)
            nu[ii] = tmp[ii] - (*coeff)*ETA_bar[i*p_dim+ii];
        
        matrixMultiply(p_dim, 1, p_dim, 1, 2, 0, nu, nu, (nu2+i));
    }    
}
