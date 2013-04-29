#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

double* rhs_z(double* rhs, double* u, double* sub_diag_z, double* hyp_diag_z, size_t N);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *u_mid, *u, *sub_diag_z, *hyp_diag_z, *diag_z, *f_cur, *f_next, *rhs;
	double h2;

	size_t N, N2, i, j, k;

	if( 1 != nlhs || nrhs != 8)
		mexErrMsgTxt("wrong number of arguments");

	u_mid= mxGetPr(prhs[0]);
	u= mxGetPr(prhs[1]);

	sub_diag_z= mxGetPr(prhs[2]);
	diag_z= mxGetPr(prhs[3]);
	hyp_diag_z= mxGetPr(prhs[4]);
	f_cur= mxGetPr(prhs[5]);
	f_next= mxGetPr(prhs[6]);

	N= mxGetM(prhs[0]);
	N2=N*N;
	h2= 1.0/((N-1)*(N-1));

	/* Output */
	plhs[0]= prhs[7];
	rhs= mxGetPr(plhs[0]);

	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			for(i= 0; i< N; ++i){
				rhs[i+ j*N+ k*N2]= h2*u_mid[i+ j*N+ k*N2]+ diag_z[i+ j*N+ k*N2]*u[i+ j*N+ k*N2]+
					0.5*(f_next[i+ j*N+ k*N2]-f_cur[i+ j*N+ k*N2]);
			}
		}
	}
	
	rhs= rhs_z(rhs, u, sub_diag_z, hyp_diag_z, N);

	return;
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
