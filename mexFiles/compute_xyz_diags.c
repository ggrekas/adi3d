#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

double * derivative_x(double * f, mwSize N);
double * derivative_y(double * f, mwSize N);
double * derivative_z(double * f, mwSize N);
double * derivative_xx(double * f, mwSize N);
double * derivative_yy(double * f, mwSize N);
double * derivative_zz(double * f, mwSize N);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *a, *Cg, *g, *Cphi, *phi, *C, *h; //input
	char *xyz; //input
	double *sub_diag, *diag, *hyp_diag; //output
	
	
	//internal variables
	double *a_d, *g_d, *phi_d, hI2, *u_d_coeff, hh, *g_dd, *phi_dd;
	
	mwSize N, N2, i, j, k, ndim;
	const mwSize *dims;

	if( 3 != nlhs || nrhs != 7)
        mexErrMsgTxt("not correct number of input arguments");

	//input
	a= mxGetPr(prhs[0]);
	Cg= mxGetPr(prhs[1]);
	g= mxGetPr(prhs[2]);
	Cphi= mxGetPr(prhs[3]);
	phi= mxGetPr(prhs[4]);
	C= mxGetPr(prhs[5]);
	xyz= (char *) mxGetData(prhs[6]);
	
	//output
	ndim = mxGetNumberOfDimensions(prhs[0]);
	dims = mxGetDimensions(prhs[0]);
	plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
	sub_diag= mxGetPr(plhs[0]);
	plhs[1] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
	diag= mxGetPr(plhs[1]);
	plhs[2] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
	hyp_diag= mxGetPr(plhs[2]);
	
	N= (mwSize) mxGetM(prhs[0]);
	N2= N*N;
	hI2= 0.5/(N-1);
	hh= 1.0/((N-1)*(N-1));

	u_d_coeff= (double *) malloc(N*N*N*sizeof(double)); //TODO free me!!!!!	
	if(xyz[0]=='x'){
		a_d= derivative_x(a, N);
		g_d= derivative_x(g, N);
		phi_d= derivative_x(phi, N);
		
		g_dd= derivative_xx(g,N);
		phi_dd= derivative_xx(phi,N);
	}else if(xyz[0]=='y'){
		a_d= derivative_y(a, N);
		g_d= derivative_y(g, N);
		phi_d= derivative_y(phi, N);
		
		g_dd= derivative_yy(g,N);
		phi_dd= derivative_yy(phi,N);
	}else if(xyz[0]=='z'){
		a_d= derivative_z(a, N);
		g_d= derivative_z(g, N);
		phi_d= derivative_z(phi, N);
		
		g_dd= derivative_zz(g,N);
		phi_dd= derivative_zz(phi,N);
	}else{
		mexErrMsgTxt("xyz argument has invalid value!\n");
	}
	
	for(k= 0; k< N; ++k){
			for(j= 0; j< N; ++j){
				for(i= 0; i< N; ++i){
					u_d_coeff[i+ j*N+ k*N2]= hI2*(a_d[i+ j*N+ k*N2]+
									Cg[i+ j*N+ k*N2]*g_d[i+ j*N+ k*N2]+
									Cphi[i+ j*N+ k*N2]*phi_d[i+ j*N+ k*N2]);
				}
			}
		}
		for(k= 0; k< N; ++k){
			for(j= 0; j< N; ++j){
				for(i= 0; i< N; ++i){
					hyp_diag[i+ j*N+ k*N2]= a[i+ j*N+ k*N2]+ u_d_coeff[i+ j*N+ k*N2];
				}
			}
		}
		for(k= 0; k< N; ++k){
			for(j= 0; j< N; ++j){
				for(i= 0; i< N; ++i){
					sub_diag[i+ j*N+ k*N2]= a[i+ j*N+ k*N2]- u_d_coeff[i+ j*N+ k*N2];
				}
			}
		}
		for(k= 0; k< N; ++k){
			for(j= 0; j< N; ++j){
				for(i= 0; i< N; ++i){
					diag[i+ j*N+ k*N2]= -2*a[i+ j*N+ k*N2]+ hh*(
									Cg[i+ j*N+ k*N2]*g_dd[i+ j*N+ k*N2]+
									Cphi[i+ j*N+ k*N2]*phi_dd[i+ j*N+ k*N2]+
									C[i+ j*N+ k*N2]/3.0);
				}
			}
		}

	free(a_d);
	free(g_d);
	free(phi_d);
	free(u_d_coeff);
	free(g_dd);
	free(phi_dd);
	return;
}


double * derivative_x(double * f, mwSize N){
	double *df, Ih;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih=0.5*(N-1); // 0.5/h
	
	df= (double *) calloc(N*N*N,sizeof(double)); //TODO free me!!!!!
	
	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			for(i= 1; i< N-1; ++i){
				df[i+ j*N+ k*N2]= Ih*(f[i+1+ j*N+ k*N2]- f[i-1+ j*N+ k*N2]);
			}
		}
	}
	return df;
}


double * derivative_y(double * f, mwSize N){
	double *df, Ih;
	mwSize N2, i,j,k;
		
	df= (double *) calloc(N*N*N,sizeof(double)); //TODO free me!!!!!

	N2=N*N;
	Ih= 0.5*(N-1); // 0.5/h
	
	for(k= 0; k< N; ++k){
		for(j= 1; j< N-1; ++j){
			for(i= 0; i< N; ++i){
				df[i+ j*N+ k*N2]= Ih*(f[i+ (j+1)*N+ k*N2]- + f[i+ (j-1)*N+ k*N2]);
			}
		}
	}
	return df;
}

double * derivative_z(double * f, mwSize N){
	double *df, Ih;
	mwSize N2, i,j,k;
	
	df= (double *) calloc(N*N*N,sizeof(double)); //TODO free me!!!!!

	N2=N*N;
	Ih= 0.5*(N-1); // 0.5/h
	
	for(k= 1; k< N-1; ++k){
		for(j= 0; j< N; ++j){
			for(i= 0; i< N; ++i){
				df[i+ j*N+ k*N2]= Ih*(f[i+ j*N+ (k+1)*N2]- f[i+ j*N+ (k-1)*N2]);
			}
		}
	}
	return df;
}
	
double * derivative_xx(double * f, mwSize N){
	double *df, Ih2;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih2= ((N-1)*(N-1)); // 1/h^2
	
	df= (double *) malloc(N*N*N*sizeof(double)); //TODO free me!!!!!
	
	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			for(i= 1; i< N-1; ++i){
				df[i+ j*N+ k*N2]= Ih2*(f[i+1+ j*N+ k*N2]- 2* f[i+ j*N+ k*N2]+ f[i-1+ j*N+ k*N2]);
			}
		}
	}
	
	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			df[j*N+ k*N2]= Ih2*2*(f[1+ j*N+ k*N2]- f[j*N+ k*N2]);
		}
	}
	
	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			df[N-1+ j*N+ k*N2]= Ih2*2*(f[N-2 + j*N+ k*N2]- f[N-1+ j*N+ k*N2]);
		}
	}
	return df;
}

double * derivative_yy(double * f, mwSize N){
	double *df, Ih2;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih2= ((N-1)*(N-1)); // 1/h^2
	
	df= (double *) malloc(N*N*N*sizeof(double)); //TODO free me!!!!!
	for(k= 0; k< N; ++k){
		for(j= 1; j< N-1; ++j){
			for(i= 0; i< N; ++i){
				df[i+ j*N+ k*N2]= Ih2*(f[i+ (j+1)*N+ k*N2]- 2* f[i+ j*N+ k*N2]+ f[i+ (j-1)*N+ k*N2]);
			}
		}
	}
	
	for(k= 0; k< N; ++k){
		for(i= 0; i< N; ++i){
			df[i+ k*N2]= Ih2*2*(f[i+ N+ k*N2]- f[i+ k*N2]);
		}
	}
	
	for(k= 0; k< N; ++k){
		for(i= 0; i< N; ++i){
			df[i+ (N-1)*N+ k*N2]= Ih2*2*(f[i+ (N-2)*N+ k*N2]- f[i+ (N-1)*N+ k*N2]);
		}
	}
	return df;
}

double * derivative_zz(double * f, mwSize N){
	double *df, Ih2;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih2= ((N-1)*(N-1)); // 1/h^2
	
	df= (double *) malloc(N*N*N*sizeof(double)); //TODO free me!!!!!
	for(k= 1; k< N-1; ++k){
		for(j= 0; j< N; ++j){
			for(i= 0; i< N; ++i){
				df[i+ j*N+ k*N2]= Ih2*(f[i+ j*N+ (k+1)*N2]- 2* f[i+ j*N+ k*N2]+ f[i+ j*N+ (k-1)*N2]);
			}
		}
	}
	
	for(j= 0; j< N; ++j){
		for(i= 0; i< N; ++i){
			df[i+ j*N]= Ih2*2*(f[i+ j*N+ N2]- f[i+ j*N]);
		}
	}
	
	for(j= 0; j< N; ++j){
		for(i= 0; i< N; ++i){
			df[i+ j*N+ (N-1)*N2]= Ih2*2*(f[i + j*N+ (N-2)*N2]- f[i+ j*N+ (N-1)*N2]);
		}
	}
	return df;
}

