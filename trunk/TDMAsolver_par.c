#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"
#include <omp.h>


void solveMatrix (int n, const double *a, double *b, const double *c,
        double *v, double *x);
void init_param(double* dest, const double* source, int N);
int par_mat_init(double ***d, double ***rhs, int N);
void par_mat_free(double **d, double **rhs, int nthreads);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
   double *a, *diag, *c, *V, *U;
   double *sub_diag, *hyp_diag;
   double **d, **rhs, b;
   int N, i, step, dims, j; /*can be replaced with mxGetM*/
   int tid, nthreads;

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
  
   #if 0
   if( N==mxGetN(prhs[1]) )
      ismatrix=1;
   else
      ismatrix=0;
   #endif
   nthreads = par_mat_init(&d, &rhs, N);

   dims = mxGetNumberOfDimensions( prhs[0] );
   #pragma omp parallel shared(N, d, rhs, sub_diag, diag, hyp_diag, U)\
           private(i, j, tid)
   {
	  tid = omp_get_thread_num();
     #pragma omp for schedule(guided) nowait
     for(i =0; i < N; ++i)
     {
		   for( j =0; j < N; ++j)
         { 
			   init_param( d[tid], diag+ j*N +i*N*N, N);
			   init_param( rhs[tid], V + j*N + i*N*N, N);
			   solveMatrix(N, sub_diag+ j*N + i*N*N, d[tid], hyp_diag+j*N+i*N*N,
				   rhs[tid], U + j*N + i*N*N);
         }
     }

   }

   par_mat_free(d, rhs, nthreads);
   
   return;
}

int par_mat_init(double ***d, double ***rhs, int N){
   int tid, nthreads, i;

#pragma omp parallel private(tid)
{
   /* Obtain and print thread id */
   tid = omp_get_thread_num();
   /* Only master thread does this */
   if (tid == 0){
      nthreads = omp_get_num_threads();
   }
}

*d= (double **)mxMalloc( (nthreads)*sizeof(double*));
*rhs= (double **)mxMalloc( (nthreads)*sizeof(double*));
  
for(i=0; i<nthreads; ++i){
   *(*d+i)= (double *)mxMalloc(N*sizeof(double));
   *(*rhs+i)= (double *)mxMalloc(N*sizeof(double));
}

return nthreads;
}

void par_mat_free(double **d, double **rhs, int nthreads){
   int i;

   for(i=0; i<nthreads; ++i){
      mxFree( *(d+i) );
      mxFree( *(rhs+i) );
   }
   mxFree(d);
   mxFree(rhs);
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