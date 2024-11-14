#include "params.h"

void Connections(int a, int b, double *gammaTil_ij, double *gammaTil_IJ, double *GammaTil_Ijk){

	int ab = a*N+b;
	int i,j,k,l;
	int ijk, il, jk, jl, kl;
	double term1, term2, term3;

	for(i=0; i<Dm; i++)	
		for(j=0; j<Dm; j++)	
			for(k=0; k<Dm; k++){

				jk = j*Dm+k;

				ijk = (i*Dm+j)*Dm+k;

				term1 = term2 = term3 = 0.0;
				
				for(l=0; l<Dm; l++){

					il = i*Dm+l;
					kl = k*Dm+l;	
					jl = j*Dm+l;

					term1 += gammaTil_IJ[ab*rank2+il]*Del(2, j, kl, a,b, gammaTil_ij);

					term2 += gammaTil_IJ[ab*rank2+il]*Del(2, k, jl, a,b, gammaTil_ij);

					term3 += gammaTil_IJ[ab*rank2+il]*Del(2, l, jk, a,b, gammaTil_ij);
				}

				GammaTil_Ijk[ijk] = 0.5*(term1 + term2 - term3);
			}

}

void ConformTensor(int a, int b, double *chi, double *gammaTil_ij, double *gammaTil_IJ, double *C_Ijk){

	int ab = a*N+b;
	int i, j, k, l;
	int ijk, il, jk; 
	double term1, term2, term3, d3;

	for(i=0; i<Dm; i++)
		for(j=0; j<Dm; j++)
			for(k=0; k<Dm; k++){

				jk = j*Dm+k;

				ijk = (i*Dm+j)*Dm+k;

				term1 = (i == j)? Del(0, k, 0, a,b, chi): 0.0;
				
				term2 = (i == k)? Del(0, j, 0, a,b, chi): 0.0;

				term3 = 0.0;
				for(l=0; l<Dm; l++){

					il = i*Dm+l;

					d3 = Del(0, l, 0, a,b, chi);
					
					term3 += gammaTil_ij[ab*rank2+jk]*gammaTil_IJ[ab*rank2+il]*d3;
				}						

				C_Ijk[ijk] = -(0.5/chi[ab])*(term1 + term2 - term3); 


				if(isinf(C_Ijk[ijk])){ printf("\nC_Ijk is Inf !\n"); exit(EXIT_FAILURE); }
                                else if(isnan(C_Ijk[ijk])){ printf("\nC_Ijk is NaN !\n"); exit(EXIT_FAILURE); }

			} // for(k)

}

int epsilon(int i, int j, int k){

        int e, ip,im, jp,jm, kp,km;

        ip = i+1;  im = i-1;
        jp = j+1;  jm = j-1;
        kp = k+1;  km = k-1;

        if(ip > 2) ip = 0;
        else if(im < 0) im = 2;

        if(jp > 2) jp = 0;
        else if(jm < 0) jm = 2;

        if(kp > 2) kp = 0;
        else if(km < 0) km = 2;

        if(i==kp && j==ip && k==jp) e = 1;
        else if((i==km) && (j==im) && (k==jm)) e = -1;
        else e = 0;

        return e;
}

