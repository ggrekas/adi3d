#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

double* rhs_x(double *rhs,double * u,double * sub_diag_x,double * hyp_diag_x, size_t N);
double* rhs_y(double* rhs, double* u, double* sub_diag_y, double* hyp_diag_y, size_t N);
double* rhs_z(double* rhs, double* u, double* sub_diag_z, double* hyp_diag_z, size_t N);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *u, *f, *sub_diag_x, *hyp_diag_x, *diag_x, *sub_diag_y, *hyp_diag_y, *diag_y,
					*sub_diag_z, *hyp_diag_z, *diag_z; //input
	double *rhs; //output
	
	double h2;

	size_t N, N2, i, j, k;
	mwSize ndim;
	const mwSize *dims;

	if( 1 != nlhs || nrhs != 11)
		mexErrMsgTxt("wrong number of arguments");

	u= mxGetPr(prhs[0]);
	f= mxGetPr(prhs[1]);

	sub_diag_x= mxGetPr(prhs[2]);
	diag_x= mxGetPr(prhs[3]);
	hyp_diag_x= mxGetPr(prhs[4]);

	sub_diag_y= mxGetPr(prhs[5]);
	diag_y= mxGetPr(prhs[6]);
	hyp_diag_y= mxGetPr(prhs[7]);

	sub_diag_z= mxGetPr(prhs[8]);
	diag_z= mxGetPr(prhs[9]);
	hyp_diag_z= mxGetPr(prhs[10]);

	N= mxGetM(prhs[0]);
	N2=N*N;
	h2= 1.0/((N-1)*(N-1)); // h^2

	//output
	ndim = mxGetNumberOfDimensions(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
	rhs= mxGetPr(plhs[0]);

	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			for(i= 0; i< N; ++i){
				rhs[i+ j*N+ k*N2]= (h2+diag_x[i+ j*N+ k*N2]+diag_y[i+ j*N+ k*N2]+diag_z[i+ j*N+ k*N2])*
					u[i+ j*N+ k*N2]+ f[i+ j*N+ k*N2];
			}
		}
	}

	rhs= rhs_x(rhs, u, sub_diag_x, hyp_diag_x, N);
	rhs= rhs_y(rhs, u, sub_diag_y, hyp_diag_y, N);
	rhs= rhs_z(rhs, u, sub_diag_z, hyp_diag_z, N);

	return;
}

double* rhs_x(double *rhs,double * u,double * sub_diag_x,double * hyp_diag_x, size_t N){
	size_t i,j,k, N2=N*N;

	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			for(i= 1; i< N-1; ++i){
				rhs[i+ j*N+ k*N2]+= sub_diag_x[i+ j*N+ k*N2]*u[i-1+ j*N+ k*N2]+
					hyp_diag_x[i+ j*N+ k*N2]*u[i+1+ j*N+ k*N2];
			}
			rhs[j*N+ k*N2]+= (hyp_diag_x[j*N+ k*N2]+ sub_diag_x[j*N+ k*N2])
				*u[1+ j*N+ k*N2];
			rhs[N-1+ j*N+ k*N2]+= (hyp_diag_x[N-1+ j*N+ k*N2]+ sub_diag_x[N-1+ j*N+ k*N2])
				*u[N-2+ j*N+ k*N2];
		}
	}

	return rhs;
}

double* rhs_y(double* rhs, double* u, double* sub_diag_y, double* hyp_diag_y, size_t N){
	size_t i,j,k, N2=N*N;

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

double* rhs_z(double* rhs, double* u, double* sub_diag_z, double* hyp_diag_z, size_t N){
	size_t i,j,k, N2=N*N;

	for(k= 1; k< N-1; ++k){
		for(j= 0; j< N; ++j){
			for(i= 0; i< N; ++i){
				rhs[i+ j*N+ k*N2]+= sub_diag_z[i+ j*N+ k*N2]*u[i+ j*N+ (k-1)*N2]+
					hyp_diag_z[i+ j*N+ k*N2]*u[i+ j*N+ (k+1)*N2];
			}
		}
	}

	/*Neumann boundary conditions*/
	for(j= 0; j< N; ++j){
		for(i= 0; i< N; ++i){
			rhs[i+ j*N]+= (hyp_diag_z[i+ j*N]+ sub_diag_z[i+ j*N])
				*u[i+ j*N+ N2];
			rhs[i+ j*N+ (N-1)*N2]+= (hyp_diag_z[i+ j*N+ (N-1)*N2]+ sub_diag_z[i+ j*N+ (N-1)*N2])
				*u[i+ j*N+ (N-2)*N2];
		}
	}
	
	return rhs;
}

