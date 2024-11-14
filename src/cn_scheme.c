#include "params.h"

extern int Eq;
extern double weight, zeta;

void CN_scheme(int rank, int i, int j, int a, int b, double *alpha, double *K, double *beta_I, double *gammaTil_IJ, double *g, double *f){

	int ij = i*Dm+j;
	int ab = a*N+b;
	int w = pow(Dm,rank);
	double f_dot, g_dot, L_beta = 0.0;
	double xi_f = 0.0, xi_g = 0.0;

	//__________________________________________________________________________
	if(weight != 1) L_beta = Lie(rank, i,j, a,b, beta_I, f);

	f_dot = L_beta + zeta + func(i,j, a,b, alpha, K, beta_I, gammaTil_IJ, f);

	if(eps != 0) xi_f = KO4th(rank, ij, a,b, f);

	g[ab*w+ij] = f[ab*w+ij] + dt*(f_dot + xi_f);


	//__________________________________________________________________________
	if(weight != 1) L_beta = Lie(rank, i,j, a,b, beta_I, g);

	g_dot = L_beta + zeta + func(i,j, a,b, alpha, K, beta_I, gammaTil_IJ, g);

	if(eps != 0) xi_g = KO4th(rank, ij, a,b, g);

	g[ab*w+ij] = f[ab*w+ij] + 0.5*dt*(g_dot + xi_g + f_dot + xi_f);


	//__________________________________________________________________________
	if(weight != 1) L_beta = Lie(rank, i,j, a,b, beta_I, g);

	g_dot = L_beta + zeta + func(i,j, a,b, alpha, K, beta_I, gammaTil_IJ, g);

	if(eps != 0) xi_g = KO4th(rank, ij, a,b, g);

	f[ab*w+ij] += 0.5*dt*(g_dot + xi_g + f_dot + xi_f);  
}

double func(int i, int j, int a, int b, double *alpha, double *K, double *beta_I, double *gammaTil_IJ, double *f){

	int ij = i*Dm+j;
	int ab = a*N+b;
	int k,l, ik, jl, kl;
	double F, term = 0.0;

	if(Eq == 1){

		for(k=0; k<Dm; k++)
			term += Del(1, k, k, a,b, beta_I);
		
		F = (2/3.)*f[ab]*(alpha[ab]*K[ab] - term);
	}

	else if(Eq == 2) F = -2*K[ab]*f[ab];

	else if(Eq == 3) F = 0.0;

	else if(Eq == 4) F = (1/3.)*alpha[ab]*pow(f[ab],2);

	else if(Eq == 5) F = 0.0;

	else if(Eq == 6){

		for(k=0; k<Dm; k++)
			term += beta_Del(1, k, j, a,b, beta_I, f);

		F = term - eta*f[ab*rank1+j];
	}

	else if(Eq == 7){

		for(k=0; k<Dm; k++)
			term += beta_Del(1, k, j, a,b, beta_I, f);

		F = term;
	}

	else if(Eq == 8){

		for(k=0; k<Dm; k++){	

			ik = i*Dm+k;

			for(l=0; l<Dm; l++){

				jl = j*Dm+l;
				kl = k*Dm+l;

				term += gammaTil_IJ[ab*rank2+kl]*f[ab*rank2+ik]*f[ab*rank2+jl];
			}	
		}

		F = alpha[ab]*(K[ab]*f[ab*rank2+ij] - 2*term);
	}

	return F;
}
