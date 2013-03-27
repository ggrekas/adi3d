#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"


void solveMatrix (int n, const double *a, double *b, const double *c,
        double *v, double *x);
void init_param(double* dest, const double* source, int N);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
   double *a, *diag, *c, *V, *U;
   double *sub_diag, *hyp_diag;
   double *d, *rhs, b;
   int N, i, step, ismatrix, dims; /*can be replaced with mxGetM*/
   
   if( 0 == nlhs || nrhs < 5)
      mexErrMsgTxt("not enough input arguments");
   plhs[0] =prhs[0];
   
   U= mxGetPr(plhs[0]);
   sub_diag= mxGetPr(prhs[1]);
   diag= mxGetPr(prhs[2]);
   hyp_diag= mxGetPr(prhs[3]);
   V= mxGetPr(prhs[4]);
   
   N= mxGetM( prhs[0] );
   
   d= (double *)mxMalloc(N*sizeof(double));
   rhs= (double *)mxMalloc(N*sizeof(double));
   
   if( N==mxGetN(prhs[1]) )
      ismatrix=1;
   else
      ismatrix=0;
   
   dims = mxGetNumberOfDimensions( prhs[0] );
   for(i =0; i < pow(N, dims-1 ); ++i){
      init_param( d, diag+i*N, N);
      init_param( rhs, V+i*N, N);
      solveMatrix(N, sub_diag+i*N*ismatrix, d, hyp_diag+i*N*ismatrix, rhs, U +i*N);
   }
   mxFree(d);
   mxFree(rhs);
   
   return;
}

void init_param(double* dest, const double* source, int N){
   int i;
   for(i=0; i<N; ++i)
      dest[i]= source[i];
   return;
}


void solveMatrix (int n, const double *a, double *B, const double *c, double *V, double *X)
{
   /*
    * n - number of equations
    * a - sub-diagonal (means it is the diagonal below the main diagonal)-
    * - indexed from 1..n-1
    * b - the main diagonal
    * c - sup-diagonal (means it is the diagonal above the main diagonal)-
    * - indexed from 0..n-2
    * v - right part
    * x - the answer
    */
   int i;
   double m;
   
   
   for (i = 1; i < n; i++){
      m = a[i]/B[i-1];
      B[i] = B[i] - m * c[i - 1];
      V[i] = V[i] - m*V[i-1];
   }
   X[n-1] = V[n-1]/ B[n-1];
   
   for (i = n - 2; i >= 0; i--){
      X[i] = (V[i] - c[i] *X[i+1]) / B[i];
   }
   return;
}