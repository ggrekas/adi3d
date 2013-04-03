#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

double* rhs_y(double* rhs, double* u, double* sub_diag_y, double* hyp_diag_y, size_t N);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *u_mid, *u, *sub_diag_y, *hyp_diag_y, *diag_y, *rhs;
	double h2;
	size_t N, N2, i, j, k;
	mwSize ndim;
	const mwSize *dims;
	
	if( 1 != nlhs || nrhs != 5)
		mexErrMsgTxt("wrong number of arguments");

	u_mid= mxGetPr(prhs[0]);
	u= mxGetPr(prhs[1]);

	sub_diag_y= mxGetPr(prhs[2]);
	diag_y= mxGetPr(prhs[3]);
	hyp_diag_y= mxGetPr(prhs[4]);

	N= mxGetM(prhs[0]);
	N2=N*N;
	h2= 1.0/((N-1)*(N-1)); //h^2

	//output
	ndim = mxGetNumberOfDimensions(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
	rhs= mxGetPr(plhs[0]);

	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			for(i= 0; i< N; ++i){
				rhs[i+ j*N+ k*N2]= h2*u_mid[i+ j*N+ k*N2]+ diag_y[i+ j*N+ k*N2]*u[i+ j*N+ k*N2];
			}
		}
	}
	
	rhs= rhs_y(rhs, u, sub_diag_y, hyp_diag_y, N);

	return;
}

double* rhs_y(double* rhs, double* u, double* sub_diag_y, double* hyp_diag_y, size_t N){
	int i,j,k, N2=N*N;

	for(k= 0; k< N; ++k){
		for(j= 1; j< N-1; ++j){
			for(i= 0; i< N; ++i){
				rhs[i+ j*N+ k*N2]+= sub_diag_y[i+ j*N+ k*N2]*u[i+ (j-1)*N+ k*N2]+
					hyp_diag_y[i+ j*N+ k*N2]*u[i+ (j+1)*N+ k*N2];
			}
		}
	}

	/*Neumann boundary conditions*/
	for(k= 0; k< N; ++k){
		for(i= 0; i< N; ++i){
			rhs[i+ k*N2]+= (hyp_diag_y[i+ k*N2]+ sub_diag_y[i+ k*N2])
				*u[i+ N+ k*N2];
			rhs[i+ (N-1)*N+ k*N2]+= (hyp_diag_y[i+ (N-1)*N+ k*N2]+ sub_diag_y[i+ (N-1)*N+ k*N2])
				*u[i+ (N-2)*N+ k*N2];
		}
	}
	return rhs;
}
