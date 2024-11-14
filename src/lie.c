#include "params.h"

extern double weight;

double Lie(int rank, int i, int j, int a, int b, double *beta_I, double *f){

	int ab = a*N+b;
	int k, ik, jk, ij = i*Dm+j;
	double L_beta, term1, term2, term3, term4;

	if(rank == 0){

		term1 = term2 = 0.0;

		for(k=0; k<Dm; k++){

			term1 += beta_Del(0, k, 0, a,b, beta_I, f);
	
			if(weight != 0.0)
			term2 += weight*Del(1, k, k, a,b, beta_I);
		}

		L_beta = term1 + term2;

		return L_beta;
	}
	
	else if(rank == 1){

		term1 = term2 = term3 = 0.0;

		for(k=0; k<Dm; k++){

			term1 += beta_Del(1, k, j, a,b, beta_I, f);
			
			term2 += f[ab*rank1+k]*Del(1, k, j, a,b, beta_I);

			if(weight != 0.0)
			term3 += weight*f[ab*rank1+j]*Del(1, k, k, a,b, beta_I);
		}

		L_beta = term1 - term2 + term3;

		return L_beta;
	}

	else if(rank == 2){

		term1 = term2 = term3 = term4 = 0.0;

		for(k=0; k<Dm; k++){

			ik = i*Dm+k;
			jk = j*Dm+k;

			term1 += beta_Del(2, k, ij, a,b, beta_I, f);
			
			term2 += f[ab*rank2+ik]*Del(1, j, k, a,b, beta_I);

			term3 += f[ab*rank2+jk]*Del(1, i, k, a,b, beta_I);

			if(weight != 0.0)
			term4 += weight*f[ab*rank2+ij]*Del(1, k, k, a,b, beta_I);
		}

		L_beta = term1 + term2 + term3 + term4;

		return L_beta;
	}
}
