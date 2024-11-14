#include "params.h"

int Eq; // auxiliary global variables
double weight, zeta, chi_p;

int main(){

	int t, a,b, ab;  // spacetime indexes
	int i,j,k,l; // tensor indexes

	double R, RTil; // (3)^(Ricci scalar)
	double *chi, *alpha, *K;
	double *beta_I, *GammaTil_I;
	double *gamma_ij, *gamma_IJ, *gammaTil_ij, *gammaTil_IJ;
	double *ATil_ij, *ATil_IJ, *R_ij;
	double *GammaTil_Ijk, *C_Ijk;	

//      conformal factor  |  lapse function  |  mean curvature
	chi = Al(rank0Dm);  alpha = Al(rank0Dm);  K = Al(rank0Dm); 

//	shift function  |  conf. contracted connection
	beta_I = Al(rank1Dm);  GammaTil_I = Al(rank1Dm);
	
//	spatial metric  |  inverse of gamma_ij  |  conformal spatial metric  |  inverse of gammaTil_ij 
	gamma_ij = Al(rank2Dm);  gamma_IJ = Al(rank2Dm);  gammaTil_ij = Al(rank2Dm);  gammaTil_IJ = Al(rank2Dm);

//	extrinsic curvature's traceless part  |  inverse of ATil_ij  |  (3)^(Ricci tensor)
	ATil_ij = Al(rank2Dm);  ATil_IJ = Al(rank2);  R_ij = Al(rank2);

//	conformal connections  |  conformal tensor
	GammaTil_Ijk = Al(rank3);  C_Ijk = Al(rank3);

// auxiliary variables _________________________________________________________________________________________________________________

	int count = 1; // interaction counter
	int ij, ik, jk, jl, kl, ijk, ikl, kij;
	
	double d1, d2;
	double term_a, term_b;
	double term1, term2, term3, term4, term5, term6;
	
	double AA, tr; 
	double Rchi_ij, RTil_ij;
	
	double b_I[Dm], DDa[rank2], A_tf[rank2];

	double *x_I = Al(Nt*Dm*2);
	double *B_I = Al(rank1Dm);
	double *g1, *g2, *g3, *g4, // chi, alpha, gammaTil_ij, K
	       *g5, *g6, *g7, *g8; // GammaTil_I, B_I, beta_I, ATil_ij

	g1 = Al(rank0Dm);  g2 = Al(rank0Dm);  g3 = Al(rank2Dm);  g4 = Al(rank0Dm);
	g5 = Al(rank1Dm);  g6 = Al(rank1Dm);  g7 = Al(rank1Dm);  g8 = Al(rank2Dm);

//_____ initial conditions _____________________________________________________________________________________________________________
	Ic(chi, alpha, K, x_I, beta_I, GammaTil_I, gamma_ij, gamma_IJ, gammaTil_ij, gammaTil_IJ, ATil_ij); 

	for(a=0; a<N; a++)
		for(b=0; b<N; b++){

			ab = a*N+b;

			g1[ab] = chi[ab];
			g2[ab] = alpha[ab];
			g4[ab] = K[ab];

			for(i=0; i<Dm; i++){

				g5[ab*rank1+i] = GammaTil_I[ab*rank1+i];
				g6[ab*rank1+i] = B_I[ab*rank1+i];
				g7[ab*rank1+i] = beta_I[ab*rank1+i];

				for(j=0; j<Dm; j++){

					ij = i*Dm+j;

					g3[ab*rank2+ij] = gammaTil_ij[ab*rank2+ij];
					g8[ab*rank2+ij] = ATil_ij[ab*rank2+ij];
				}
			}
		}
	

	if(key_dynamics == 1) //________________________________________________________________________________________________________
	for(t=0; t<Nt; t++){  

		for(a=0; a<N; a++)
			for(b=0; b<N; b++){

				ab = a*N+b;
	
			// chi: ·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=

				Eq = 1;

				weight = 0.0;

				zeta = 0.0;

				CN_scheme(0, 0,0, a,b, alpha, K, beta_I, gammaTil_IJ, g1, chi);

			//	Diverge(t, ab, chi, "chi");
	
			// alpha: 1+log slicing ·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=

				Eq = 2;

				weight = 0.0;

				zeta = 0.0;

				CN_scheme(0, 0,0, a,b, alpha, K, beta_I, gammaTil_IJ, g2, alpha);
	
			//	Diverge(t, ab, alpha, "alpha");

			// gammaTil_ij: ·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=

				Eq = 3;

				weight = -2/3.;

				for(i=0; i<Dm; i++)
					for(j=0; j<Dm; j++){
					
						ij = i*Dm+j;
			
						zeta = -2*alpha[ab]*ATil_ij[ab*rank2+ij];

						CN_scheme(2, i,j, a,b, alpha, K, beta_I, gammaTil_IJ, g3, gammaTil_ij);
	
			//			Diverge(t, ab*rank2+ij, gammaTil_ij, "gammaTil_ij");
					}

				// _____________________________________________________________________________________________________
				Metric(t, a,b, chi, x_I, beta_I, gamma_ij, gamma_IJ, gammaTil_ij, gammaTil_IJ);

				Connections(a,b, gammaTil_ij,gammaTil_IJ, GammaTil_Ijk);

				ConformTensor(a,b, chi, gammaTil_ij,gammaTil_IJ, C_Ijk);

				// ATil_IJ _____________________________________________________________________________________________
				for(i=0; i<Dm; i++)
					for(j=0; j<Dm; j++){

						ij = i*Dm+j;

						ATil_IJ[ij] = 0.0;							
						for(k=0; k<Dm; k++){

							ik = i*Dm+k;	
	
							for(l=0; l<Dm; l++){

								jl = j*Dm+l;
								kl = k*Dm+l;

								ATil_IJ[ij] += gammaTil_IJ[ab*rank2+ik]*gammaTil_IJ[ab*rank2+jl]
									     * ATil_ij[ab*rank2+kl];
							}
						}
					}

			// K: ·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=

				Eq = 4;

				weight = 0.0;

				term_a = Del2(0, 0, a,b, alpha);

				term_b = 0.0;
				for(i=0; i<Dm; i++){

					d1 = Del(0, i, 0, a,b, chi);

					for(j=0; j<Dm; j++){
		
						ij = i*Dm+j;
					
						d2 = Del(0, j, 0, a,b, alpha);

						term_b += (1./chi[ab])*d1*gamma_IJ[ab*rank2+ij]*d2;
					}
				}

				term1 = term_a - 1.5*term_b;  

				AA = 0.0;
				for(i=0; i<Dm; i++)
					for(j=0; j<Dm; j++){
	
						ij = i*Dm+j;

						AA += ATil_ij[ab*rank2+ij]*ATil_IJ[ij];
					}

	
				zeta = -term1 + alpha[ab]*AA;

				CN_scheme(0, 0,0, a,b, alpha, K, beta_I, gammaTil_IJ, g4, K);

			//	Diverge(t, ab, K, "K");	

			// GammaTil_I:·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=

				Eq = 5;

				weight = 2/3.;

				for(i=0; i<Dm; i++){

					term1 = term2 = term3 = term4 = term5 = 0.0;
			
					for(j=0; j<Dm; j++){			

						ij = i*Dm+j;

						term1 += ATil_IJ[ij]*Del(0, j, 0, a,b, alpha);

						term2 += (1./chi[ab])*ATil_IJ[ij]*Del(0, j, 0, a,b, chi);

						term3 += gammaTil_IJ[ab*rank2+ij]*Del(0, j, 0, a,b, K);
					
						for(k=0; k<Dm; k++){
											
							jk = j*Dm+k; 

							ijk = ij*Dm+k; 

							term4 += GammaTil_Ijk[ijk]*ATil_IJ[jk]; 

							term5 += gammaTil_IJ[ab*rank2+ij]*Ddel(1, j,k, k, a,b, beta_I);


							Diverge(t, ijk, GammaTil_Ijk, "GammaTil_Ijk"); 
						}
					}
				
					term6 = (1./chi[ab])*Del2(1, i, a,b, beta_I);

					zeta = -2*term1 + 2*alpha[ab]*(term4 - 1.5*term2 - (2/3.)*term3) + (1/3.)*term5 + term6;

					b_I[i] = zeta; // auxiliary: Gamma-driver	

					CN_scheme(1, 0,i, a,b, alpha, K, beta_I, gammaTil_IJ, g5, GammaTil_I);
	
			//		Diverge(t, ab*rank1+i, GammaTil_I, "GammaTil_I");
				}

			// beta_I: Gamma-driver ·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=
	
				for(i=0; i<Dm; i++){
			
					Eq = 6;

					weight = 2/3.;
					term1 = Lie(1, 0,i, a,b, beta_I, GammaTil_I); 

					term2 = 0.0;
					for(j=0; j<Dm; j++) 
						term2 += beta_Del(1, j, i, a,b, beta_I, GammaTil_I);
	
					weight = 1;

					zeta = b_I[i] + term1 - term2;
	
					CN_scheme(1, 0,i, a,b, alpha, K, beta_I, gammaTil_IJ, g6, B_I);

			//		Diverge(t, ab*rank1+i, B_I, "B_I");

				} // for(i)
				// _____________________________________________________________________________________________________
				for(i=0; i<Dm; i++){

					Eq = 7;

					weight = 1; 
			
					zeta = (3/4.)*B_I[ab*rank1+i];

					CN_scheme(1, 0,i, a,b, alpha, K, beta_I, gammaTil_IJ, g7, beta_I);

			//		Diverge(t, ab*rank1+i, beta_I, "beta_I");
				}

			// R_ij: =·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=

				R = AA - (2/3.)*pow(K[ab],2);
					
				for(i=0; i<Dm; i++)
					for(j=0; j<Dm; j++){

						ij = i*Dm+j;

						RTil_ij = Ricci_til(i,j, a,b, chi, gamma_ij, gamma_IJ, GammaTil_I, GammaTil_Ijk);

						Rchi_ij = Ricci_chi(i,j, a,b, chi, gammaTil_ij, gammaTil_IJ, GammaTil_I, GammaTil_Ijk);

						R_ij[ij] = RTil_ij + Rchi_ij;

			//			Diverge(t, ij, R_ij, "R_ij");
					}
			
			// ATil_ij: ·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=·=

				Eq = 8;

				weight = -2/3.;

				tr = 0.0;
				for(i=0; i<Dm; i++)
					for(j=0; j<Dm; j++){

						ij = i*Dm+j;

						term_a = Ddel(0, i,j, 0, a,b, alpha);

						term_b = 0.0;
						for(k=0; k<Dm; k++){

							kij = (k*Dm+i)*Dm+j;

							d1 = Del(0, k, 0, a,b, alpha);

							term_b += (GammaTil_Ijk[kij] + C_Ijk[kij])*d1;
						}

						DDa[ij] = term_a - term_b; 
		
						tr += (i==j)? DDa[ij]: 0.0; // trace
					}


				for(i=0; i<Dm; i++)
					for(j=0; j<Dm; j++){

						ij = i*Dm+j;
						
						term1 = DDa[ij] - (1/3.)*tr*gamma_ij[ab*rank2+ij];

						term2 = alpha[ab]*(R_ij[ij] - (1/3.)*R*gamma_ij[ab*rank2+ij]);		
		
						zeta = chi[ab]*(term2 - term1);
		
						CN_scheme(2, i,j, a,b, alpha, K, beta_I, gammaTil_IJ, g8, ATil_ij);

			//			Diverge(t, ab*rank2+ij, ATil_ij, "ATil_ij");
					} 

			} // for(b)


		if((t % (int)(Nt/Nf)) == 0.0){

			printf("slice %d\n", count);

			OutScalar(count, chi, alpha, K);
			OutVector(count, beta_I, GammaTil_I, x_I);
			OutTensor(count, gammaTil_ij, ATil_ij);

			count++;
		}

	} // for(t)
	
//______________________________________________________________________________________________________________________________________

	free(chi);  free(alpha);  free(K);
	free(beta_I);  free(B_I);  free(GammaTil_I);  free(x_I);  
	free(gamma_ij);  free(gamma_IJ);  free(gammaTil_ij);  free(gammaTil_IJ);
	free(ATil_ij);  free(ATil_IJ);  free(R_ij);	
	free(GammaTil_Ijk);  free(C_Ijk);  

	free(g1);  free(g2);  free(g3);  free(g4);
	free(g5);  free(g6);  free(g7);  free(g8);
}
