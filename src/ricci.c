#include "params.h"

double Ricci_til(int i, int j, int a, int b, double *chi, double *gammaTil_ij, double *gammaTil_IJ, double *GammaTil_I, double *GammaTil_Ijk){

	int ab = a*N+b;
	int k, l, m, n;
	int ij,ik,il,in, jk,jl,jn, kl, mn;
	int lik,ljk, mik,mil,mjk, njk,nlm;

	double d1, d2;
	double term_a, term_b, term_c;
	double term1, term2, term3, term4;

	ij = i*Dm+j;
	
	// term1: laplacian of gammaTil_ij __________________________________________________________
	term1 = (1./chi[ab])*Del2(2, ij, a,b, gammaTil_ij);

	// term2: ___________________________________________________________________________________
	term2 = 0.0;

	for(k=0; k<Dm; k++){

		ik = i*Dm+k;
		jk = j*Dm+k;
					
		d1 = Del(1, j, k, a,b, GammaTil_I);

		d2 = Del(1, i, k, a,b, GammaTil_I);
		
		term2 += gammaTil_ij[ab*rank2+ik]*d1 + gammaTil_ij[ab*rank2+jk]*d2; 	
	}

	// term3: ___________________________________________________________________________________
	term3 = 0.0;

	for(k=0; k<Dm; k++)
		for(l=0; l<Dm; l++){

			il = i*Dm+l;
			jl = j*Dm+l;
			kl = k*Dm+l;

			ljk = (l*Dm+j)*Dm+k;
			lik = (l*Dm+i)*Dm+k;

			d1 = -Del(2, l, kl, a,b, gammaTil_IJ);

			term_a = gammaTil_ij[ab*rank2+il]*GammaTil_Ijk[ljk];

			term_b = gammaTil_ij[ab*rank2+jl]*GammaTil_Ijk[lik];

			term3 += d1*(term_a + term_b);
		}


	// term4: ___________________________________________________________________________________
	term4 = 0.0;

	for(k=0; k<Dm; k++)
		for(l=0; l<Dm; l++){

			kl = k*Dm+l;

			term_a = term_b = term_c = 0.0;

			for(m=0; m<Dm; m++){
			
				mik = (m*Dm+i)*Dm+k;
				mjk = (m*Dm+j)*Dm+k;
				mil = (m*Dm+i)*Dm+l;
				
				for(n=0; n<Dm; n++){

					in = i*Dm+n;
					jn = j*Dm+n;
					mn = m*Dm+n;
					
					njk = (n*Dm+j)*Dm+k;
					nlm = (n*Dm+l)*Dm+m;			
	
					term_a += GammaTil_Ijk[mik]*gammaTil_ij[ab*rank2+jn]*
						  GammaTil_Ijk[nlm];

					term_b += GammaTil_Ijk[mjk]*gammaTil_ij[ab*rank2+in]*
						  GammaTil_Ijk[nlm];

					term_c += GammaTil_Ijk[mil]*gammaTil_ij[ab*rank2+mn]*
						  GammaTil_Ijk[njk];
				}

			} // for(m)

			term4 += gammaTil_IJ[ab*rank2+kl]*(term_a + term_b + term_c);
		}

	return -0.5*term1 + 0.5*term2 + 0.5*term3 + term4;
}

double Ricci_chi(int i, int j, int a, int b, double *chi, double *gamma_ij, double *gamma_IJ, double *GammaTil_I, double *GammaTil_Ijk){

	int k, l;
	int ab = a*N+b;
	int ij, kl, kij;

	double d1, d2, d3;
	double term_a, term_b;
	double term1, term2, term3, term4; 

	ij = i*Dm+j;

	d1 = Del(0, i, 0, a,b, chi);
	d2 = Del(0, j, 0, a,b, chi);

	// term1: _____________________________________________	
	term_a = Ddel(0, i,j, 0, a,b, chi);

	term_b = 0.0;
	for(k=0; k<Dm; k++){

		kij = (k*Dm+i)*Dm+j;

		d3 = Del(0, k, 0, a,b, chi);

		term_b += GammaTil_Ijk[kij]*d3;
	}

	term1 = (1./chi[ab])*(term_a - term_b); 

	// term2: _____________________________________________	
	term_a = (1./chi[ab])*Del2(0, 0, a,b, chi); 

	term_b = pow(chi[ab],-2)*d1*gamma_IJ[ab*rank2+ij]*d2;

	term2 = gamma_ij[ab*rank2+ij]*(term_a - term_b);

	// term3: _____________________________________________	
	term3 = pow(chi[ab],-2)*d1*d2;

	// term4: _____________________________________________	
	term_a = 0.0;
	for(k=0; k<Dm; k++){

		d1 = Del(0, k, 0, a,b, chi);

		for(l=0; l<Dm; l++){

			kl = k*Dm+l;

			d2 = Del(0, l, 0, a,b, chi);

			term_a += d1*gamma_IJ[ab*rank2+kl]*d2;
		}
	}
	
	term4 = pow(chi[ab],-2)*gamma_ij[ab*rank2+ij]*term_a;


	return 0.5*term1 + 0.5*term2 - 0.25*term3 - 1.5*term4;
}
