#include "params.h"

void Ic(double *chi, double *alpha, double *K, double *x_I, double *beta_I, double *GammaTil_I, 
	double *gamma_ij, double *gamma_IJ, double *gammaTil_ij, double *gammaTil_IJ, double *ATil_ij){

        if(r_s <= dx){

                printf("\nSchwarzschild radius is smaller than the grid spacing!\n");
                exit(EXIT_FAILURE);
        }

        printf("\n===============================\n");
        printf(" Binary Black Hole: \n");
        printf("  →  grid's length = %gkm.\n", 1.*L);
        printf("  →  Simulation over %gs.\n", 1.*T);
        printf("  →  %gs for each snapshot.\n", 1.*T/Nf);
        printf("  →  r_s = %gkm.\n", r_s);
        printf("===============================\n");

	// ____________________________________________________________________________________________________________________________________
	int i,j,k,w, ij;
	int a,b, ab, eta_ij;
	int x,y;

	double x1,y1, r1, x2,y2, r2;

	double p, Ap, As, a_ij, AA, psi;
	double k_sqr, Norm = pow(0.01,Dm);

	double *psi_BL = Al(rank0Dm);

	double *r_b = Al(2*rank0Dm);
	double *q_b = Al(2*rank1Dm);

	double *u = Al(rank0Dm);
	double *s = Al(rank0Dm);

	fftw_complex *s_c = fftw_alloc_complex(rank0Dm);
	fftw_complex *u_c = fftw_alloc_complex(rank0Dm);

	fftw_plan s_r2c = fftw_plan_dft_r2c_1d(rank0Dm, s, s_c, FFTW_ESTIMATE);
	fftw_plan u_c2r = fftw_plan_dft_c2r_1d(rank0Dm, u_c, u, FFTW_ESTIMATE);

	// coordinates:________________________________________________________________________________________________________________________
	for(a=0; a<N; a++)
		for(b=0; b<N; b++){

			ab = a*N+b; 
		// _____________________________________
			x1 = a*dx-0.5*L;
			y1 = b*dx-0.5*(L-d*r_s);

			r1 = sqrt(pow(x1,2)+pow(y1,2));	

			r_b[ab*2+0] = r1;

			q_b[(ab*rank1+0)*2+0] = x1/r1;  
			q_b[(ab*rank1+1)*2+0] = y1/r1;

			if(r1 == 0){

				x_I[(0*Dm+0)*2+0] = a*dx;
				x_I[(0*Dm+1)*2+0] = b*dx;
			}
		// _____________________________________
			x2 = a*dx-0.5*L;
			y2 = b*dx-0.5*(L+d*r_s);
		
			r2 = sqrt(pow(x2,2)+pow(y2,2));	

			r_b[ab*2+1] = r2;

			q_b[(ab*rank1+0)*2+1] = x2/r2;
			q_b[(ab*rank1+1)*2+1] = y2/r2;

			if(r2 == 0){

				x_I[(0*Dm+0)*2+1] = a*dx;
				x_I[(0*Dm+1)*2+1] = b*dx;
			} 
		}

	
	__momentum_spin__  // Brandt-Brügmann initial data:____________________________________________________________________________________
	p = 0.5*( sqrt(pow(P[0*2+0],2)+pow(P[0*2+1],2))
		+ sqrt(pow(P[1*2+0],2)+pow(P[1*2+1],2)) );

	if(p != 0)
	for(a=0; a<N; a++)
		for(b=0; b<N; b++){

			ab = a*N+b;

			for(w=0; w<2; w++) // black-hole loop
				for(i=0; i<Dm; i++)
					for(j=0; j<Dm; j++){

						ij = i*Dm+j;

						eta_ij = (i==j)? 1: 0;

						As = (2/r_b[ab*2+w])*(S[w*2+0]*q_b[(ab*rank1+1)*2+w] - S[w*2+1]*q_b[(ab*rank1+0)*2+w]) \
						   * (q_b[(ab*rank1+i)*2+w] + q_b[(ab*rank1+j)*2+w]);

						for(k=0; k<Dm; k++){

							Ap = P[w*2+i]*q_b[(ab*rank1+j)*2+w] + P[w*2+j]*q_b[(ab*rank1+i)*2+w]
							   - (eta_ij-q_b[(ab*rank1+i)*2+w]*q_b[(ab*rank1+j)*2+w])		
							   * P[w*2+k]*q_b[(ab*rank1+k)*2+w];
			
							a_ij = 1.5*pow(r_b[ab*2+w],-2)*(Ap + As);

							if(isnan(a_ij) != 1) ATil_ij[ab*rank2+ij] += a_ij;
						}
					} 
		} // for(b)

	// u:__________________________________________________________________________________________________________________________________
	if(p != 0)
	for(a=0; a<N; a++)
		for(b=0; b<N; b++){

			ab = a*N+b;

			u[ab] = 1 + 6e-5*pow(L,2)*pow(d,0.25)*(exp(p)-1)*(1/r_b[ab*2+0] + 1/r_b[ab*2+1]);
		}
	
	// s:__________________________________________________________________________________________________________________________________
	for(a=0; a<N; a++)
		for(b=0; b<N; b++){

			ab = a*N+b; 

			psi_BL[ab] = 1 + 0.25*r_s*(1/r_b[ab*2+0] + 1/r_b[ab*2+1]);
		
			if(p != 0){	

				AA = 0.0;
				for(i=0; i<Dm; i++)
					for(j=0; j<Dm; j++){

						ij = i*Dm+j;
		
						AA += (i==j)? pow(ATil_ij[ab*rank2+ij],2): 0.0;
					}
			

				s[ab] = -0.125*AA*pow(psi_BL[ab] + u[ab],-7);
			}
		}


	if(p != 0) fftw_execute(s_r2c);

	// u_c:________________________________________________________________________________________________________________________________
	if(p != 0){

		for(a=0; a<N; a++)
			for(b=0; b<N; b++){

				ab = a*N+b; 

				k_sqr = 4*pow(pi/N,2)*(a*a + b*b); 
		
				if(k_sqr == 0.0) u_c[ab] = 0.0;
				else u_c[ab] = -Norm*(1/k_sqr)*s_c[ab];	
			}


		fftw_execute(u_c2r);
	}

	// chi:________________________________________________________________________________________________________________________________
	for(a=0; a<N; a++)
		for(b=0; b<N; b++){

			ab = a*N+b; // ab = {0,1,2, ..., N^Dm}

			chi[ab] = pow(psi_BL[ab] + u[ab],-4);

			for(i=0; i<Dm; i++)	
				for(j=0; j<Dm; j++){
			
					ij = i*Dm+j; // Dm = 2 →  ij = {0·2+0, 0·2+1, 1·2+0, 1·2+1}

					//           gammaTil_ij = eta_ij
					gammaTil_ij[ab*rank2+ij] = (i==j)? 1.0: 0.0; 

					ATil_ij[ab*rank2+ij] *= pow(chi[ab],1.5);

				} // each grid point has a set of metric coefficients
				  // Dm = 2 →  ab*rank2+ij = {(0·4+0, 0·4+1, 0·4+2, 0·4+3), (1·4+0, 1·4+1, 1·4+2, 1·4+3), ..., rank2·N^Dm}
		}


	// metric:_____________________________________________________________________________________________________________________________
	for(a=0; a<N; a++)
		for(b=0; b<N; b++)
			Metric(0, a,b, chi, x_I, beta_I, gamma_ij, gamma_IJ, gammaTil_ij, gammaTil_IJ);


	// alpha:______________________________________________________________________________________________________________________________
	for(a=0; a<N; a++)
		for(b=0; b<N; b++){
		
			ab = a*N+b;
			
			if(lapse_Ic == 0) alpha[ab] = 1.0;
			else if(lapse_Ic == 1) alpha[ab] = 2*pow(chi[ab], 0.25)-1;  // 2/psi-1
			else if(lapse_Ic == 2) alpha[ab] = pow(chi[ab],0.5); 	    // psi^-2
		}


	// Output:_____________________________________________________________________________________________________________________________
	OutScalar(0, chi, alpha, K);
	OutVector(0, beta_I, GammaTil_I, x_I);
	OutTensor(0, gammaTil_ij, ATil_ij);

	// ____________________________________________________________________________________________________________________________________
	fftw_destroy_plan(s_r2c);
	fftw_destroy_plan(u_c2r);

	fftw_free(u_c);  fftw_free(s_c);	

	free(s);  free(u);

	free(psi_BL);

	free(r_b);  free(q_b);
}
