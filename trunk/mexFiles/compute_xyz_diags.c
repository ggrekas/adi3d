/*
// How to fill the tmp struct : i.e. demaloc
//======================================
// demaloc.udCoef= zeros(dim,dim,dim);
// demaloc.ad= zeros(dim,dim,dim);
// demaloc.gd= zeros(dim,dim,dim);
// demaloc.phid= zeros(dim,dim,dim);
// demaloc.gdd= zeros(dim,dim,dim);
// demaloc.phidd= zeros(dim,dim,dim);
//
// [sub_diag2, diag2, hyp_diag2]= C_compute_xyz_diags(a, Cg, g, Cphi, phi,...
// 	C, xyz, demaloc);
//====================================== 

// WARNING!!! Struct can have any name but must have the SAME field names and 
// the arrays must NOT be shared 
// ---------------THIS IS WRONG---------------
// demaloc.phidd= zeros(dim,dim,dim);
// demaloc.ad= demaloc.phidd;
// this does not create a new array, but the array is shared by the 2 pointers
// ---------------THIS IS WRONG---------------
*/
#include<math.h>
#include<mex.h>
#include<stdlib.h>
#include "matrix.h"

double * derivative_x(double * f, mwSize N, double *df);
double * derivative_y(double * f, mwSize N, double *df);
double * derivative_z(double * f, mwSize N, double *df);
double * derivative_xx(double * f, mwSize N, double *df);
double * derivative_yy(double * f, mwSize N, double *df);
double * derivative_zz(double * f, mwSize N, double *df);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *a, *Cg, *g, *Cphi, *phi, *C; /*input*/
	char *xyz; /*input*/
	double *sub_diag, *diag, *hyp_diag; /*output*/
		
	/*internal variables*/
	double *a_d, *g_d, *phi_d, hI2, *u_d_coeff, hh, *g_dd, *phi_dd;
	double phiTerm, gTerm;
	size_t dimCphi, dimCg;
	mwSize N, N2, i, j, k, ndim;
	const mwSize *dims;
	mxArray *tmp;

	if( 3 != nlhs || nrhs != 8)
		mexErrMsgTxt("not correct number of input arguments");

	/*input*/
	a= mxGetPr(prhs[0]);
	Cg= mxGetPr(prhs[1]);
	g= mxGetPr(prhs[2]);
	Cphi= mxGetPr(prhs[3]);
	phi= mxGetPr(prhs[4]);
	C= mxGetPr(prhs[5]);
	xyz= (char *) mxGetData(prhs[6]);
	
	dimCphi= mxGetM(prhs[3]);
	dimCg= mxGetM(prhs[1]);
	
	tmp = mxGetField(prhs[7],0,"udCoef");
	if( tmp == NULL )
        mexErrMsgTxt("udCoef not found");
	u_d_coeff= mxGetPr(tmp);
	
	tmp = mxGetField(prhs[7],0,"ad");
	if( tmp == NULL )
        mexErrMsgTxt("ad not found");
	a_d= mxGetPr(tmp);
	
	tmp = mxGetField(prhs[7],0,"gd");
	if( tmp == NULL )
        mexErrMsgTxt("gd not found");
	g_d= mxGetPr(tmp);
	
	tmp = mxGetField(prhs[7],0,"phid");
	if( tmp == NULL )
        mexErrMsgTxt("phid not found");
	phi_d= mxGetPr(tmp);

	tmp = mxGetField(prhs[7],0,"gdd");
	if( tmp == NULL )
        mexErrMsgTxt("gdd not found");
	g_dd= mxGetPr(tmp);
	
	tmp = mxGetField(prhs[7],0,"phidd");
	if( tmp == NULL )
        mexErrMsgTxt("phidd not found");
	phi_dd= mxGetPr(tmp);
	
	
	/*output*/
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

	if(xyz[0]=='x'){
		a_d= derivative_x(a, N, a_d);
		g_d= derivative_x(g, N, g_d);
		phi_d= derivative_x(phi, N, phi_d);
		
		g_dd= derivative_xx(g,N, g_dd);
		phi_dd= derivative_xx(phi,N, phi_dd);
	}else if(xyz[0]=='y'){
		a_d= derivative_y(a, N, a_d);
		g_d= derivative_y(g, N, g_d);
		phi_d= derivative_y(phi, N, phi_d);
		
		g_dd= derivative_yy(g,N, g_dd);
		phi_dd= derivative_yy(phi,N, phi_dd);
	}else if(xyz[0]=='z'){
		a_d= derivative_z(a, N, a_d);
		g_d= derivative_z(g, N, g_d);
		phi_d= derivative_z(phi, N, phi_d);
		
		g_dd= derivative_zz(g,N, g_dd);
		phi_dd= derivative_zz(phi,N, phi_dd);
	}else{
		mexErrMsgTxt("xyz argument has invalid value!\n");
	}
	
	for(k= 0; k< N; ++k){
			for(j= 0; j< N; ++j){
				for(i= 0; i< N; ++i){
					if (dimCphi> 1){
						phiTerm= Cphi[i+ j*N+ k*N2]*phi_d[i+ j*N+ k*N2];
					}else{
						phiTerm= Cphi[0]*phi_d[i+ j*N+ k*N2];
					}
					
					if (dimCg> 1){
						gTerm= Cg[i+ j*N+ k*N2]*g_d[i+ j*N+ k*N2];
					}else{
						gTerm= Cg[0]*g_d[i+ j*N+ k*N2];
					}
					
					u_d_coeff[i+ j*N+ k*N2]= hI2*(a_d[i+ j*N+ k*N2]+
									gTerm+ phiTerm);
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
					if (dimCphi> 1){
						phiTerm= Cphi[i+ j*N+ k*N2]*phi_dd[i+ j*N+ k*N2];
					}else{
						phiTerm= Cphi[0]*phi_dd[i+ j*N+ k*N2];
					}
					
					if (dimCg> 1){
						gTerm= Cg[i+ j*N+ k*N2]*g_dd[i+ j*N+ k*N2];
					}else{
						gTerm= Cg[0]*g_dd[i+ j*N+ k*N2];
					}
										
					diag[i+ j*N+ k*N2]= -2*a[i+ j*N+ k*N2]+ hh*(
									gTerm+ phiTerm+ C[i+ j*N+ k*N2]/3.0);
				}
			}
		}

	return;
}


double * derivative_x(double * f, mwSize N, double *df){
	double Ih;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih=0.5*(N-1); /* 0.5/h*/
	
	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			for(i= 1; i< N-1; ++i){
				df[i+ j*N+ k*N2]= Ih*(f[i+1+ j*N+ k*N2]- f[i-1+ j*N+ k*N2]);
			}
			df[j*N+ k*N2]= 0;
			df[N-1+ j*N+ k*N2]= 0;
		}
	}
	return df;
}


double * derivative_y(double * f, mwSize N, double *df){
	double Ih;
	mwSize N2, i,j,k;
		
	N2=N*N;
	Ih= 0.5*(N-1); /* 0.5/h*/
	
	for(k= 0; k< N; ++k){
		for(j= 1; j< N-1; ++j){
			for(i= 0; i< N; ++i){
				df[i+ j*N+ k*N2]= Ih*(f[i+ (j+1)*N+ k*N2]- + f[i+ (j-1)*N+ k*N2]);
			}
		}
			for(i= 0; i< N; ++i){
				df[i+ k*N2]= 0;
				df[i+ (N-1)*N+ k*N2]= 0;
			}
	}	
	return df;
}

double * derivative_z(double * f, mwSize N, double * df){
	double Ih;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih= 0.5*(N-1); /* 0.5/h*/
	
	for(k= 1; k< N-1; ++k){
		for(j= 0; j< N; ++j){
			for(i= 0; i< N; ++i){
				df[i+ j*N+ k*N2]= Ih*(f[i+ j*N+ (k+1)*N2]- f[i+ j*N+ (k-1)*N2]);
			}
		}
	}
	
	for(j= 0; j< N; ++j){
		for(i= 0; i< N; ++i){
			df[i+ j*N]= 0;
		}
	}
	
	for(j= 0; j< N; ++j){
		for(i= 0; i< N; ++i){
			df[i+ j*N+ (N-1)*N2]= 0;
		}
	}
	return df;
}
	
double * derivative_xx(double * f, mwSize N, double *df){
	double Ih2;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih2= ((N-1)*(N-1)); /* 1/h^2*/
	
	for(k= 0; k< N; ++k){
		for(j= 0; j< N; ++j){
			for(i= 1; i< N-1; ++i){
				df[i+ j*N+ k*N2]= Ih2*(f[i+1+ j*N+ k*N2]- 2* f[i+ j*N+ k*N2]+ f[i-1+ j*N+ k*N2]);
			}
			df[j*N+ k*N2]= Ih2*2*(f[1+ j*N+ k*N2]- f[j*N+ k*N2]);
			df[N-1+ j*N+ k*N2]= Ih2*2*(f[N-2 + j*N+ k*N2]- f[N-1+ j*N+ k*N2]);
		}
	}
	return df;
}

double * derivative_yy(double * f, mwSize N, double *df){
	double Ih2;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih2= ((N-1)*(N-1)); /* 1/h^2*/
	
	for(k= 0; k< N; ++k){
		for(j= 1; j< N-1; ++j){
			for(i= 0; i< N; ++i){
				df[i+ j*N+ k*N2]= Ih2*(f[i+ (j+1)*N+ k*N2]- 2* f[i+ j*N+ k*N2]+ f[i+ (j-1)*N+ k*N2]);
			}
		}
		for(i= 0; i< N; ++i){
			df[i+ k*N2]= Ih2*2*(f[i+ N+ k*N2]- f[i+ k*N2]);
			df[i+ (N-1)*N+ k*N2]= Ih2*2*(f[i+ (N-2)*N+ k*N2]- f[i+ (N-1)*N+ k*N2]);
		}
	}
	return df;
}

double * derivative_zz(double * f, mwSize N, double *df){
	double Ih2;
	mwSize N2, i,j,k;
	
	N2=N*N;
	Ih2= ((N-1)*(N-1)); /* 1/h^2*/
	
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
			df[i+ j*N+ (N-1)*N2]= Ih2*2*(f[i + j*N+ (N-2)*N2]- f[i+ j*N+ (N-1)*N2]);
		}
	}
	
	return df;
}
