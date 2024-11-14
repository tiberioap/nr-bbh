#include "params.h"

extern double chi_p;

void Metric(int t, int a, int b, double *chi, double *x_I, double *beta_I, double *gamma_ij, 
	    double *gamma_IJ, double *gammaTil_ij, double *gammaTil_IJ){

	int i,j, ij, tau;
	int ab = a*N+b;
	int ap,am, bp,bm;
	int apb,amb, abp,abm;
	int x,y, xy_1,xy_2;

	if(chi[ab] <= 0.0){

		ap = (a+1)%N;	am = (a-1+N)%N;
		bp = (b+1)%N;	bm = (b-1+N)%N;

		apb = ap*N+b;	amb = am*N+b;
		abp = a*N+bp;	abm = a*N+bm;

		chi[ab] = delta*(chi[apb]+chi[amb]+chi[abp]+chi[abm]);

		chi_p = chi[ab];
	}	
	else{
		tau = (t+1)%Nt;

		x = (int)(x_I[(t*2+0)*2+0]/dx);
		y = (int)(x_I[(t*2+1)*2+0]/dx);

		xy_1 = x*N+y;

		x = (int)(x_I[(t*2+0)*2+1]/dx);
		y = (int)(x_I[(t*2+1)*2+1]/dx);

		xy_2 = x*N+y;

		if(ab == xy_1){

			chi[ab] = chi_p;

			for(i=0; i<Dm; i++)
				x_I[(tau*Dm+i)*2+0] = x_I[(t*Dm+i)*2+0]-dt*beta_I[ab*rank1+i];

		}
		else if(ab == xy_2){

			chi[ab] = chi_p;

			for(i=0; i<Dm; i++)
				x_I[(tau*Dm+i)*2+1] = x_I[(t*Dm+i)*2+1]-dt*beta_I[ab*rank1+i];

		}
	}

	for(i=0; i<Dm; i++)
		for(j=0; j<Dm; j++){

			ij = i*Dm+j;

			gamma_ij[ab*rank2+ij] = (1./chi[ab])*gammaTil_ij[ab*rank2+ij];
		}


	Inv(t, ab, gammaTil_ij, gammaTil_IJ);
	Inv(t, ab,  gamma_ij, gamma_IJ);
}
