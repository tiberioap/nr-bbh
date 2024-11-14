#include "params.h"

void OutScalar(int t, double *chi, double *alpha, double *K){

	int a,b;
	char name[50];

	// chi:______________________________________________________________________________________
	if(out_scalar==0.1 || out_scalar==0.02 || out_scalar==0.03 || out_scalar==0.006){

		sprintf(name,"../database/chi%d.dat", t);
		FILE *scalar = fopen(name,"w");

		for(a=0; a<N; a++){
			for(b=0; b<N; b++)
				fprintf(scalar, "%g ", chi[a*N+b]);
			
			fprintf(scalar, "\n");
		}
		
		fclose(scalar);
	}
	// alpha:____________________________________________________________________________________
	if(out_scalar==0.2 || out_scalar==0.02 || out_scalar==0.06 || out_scalar==0.006){

		sprintf(name,"../database/alpha%d.dat", t);
		FILE *scalar = fopen(name,"w");

		for(a=0; a<N; a++){
			for(b=0; b<N; b++)
				fprintf(scalar, "%g ", alpha[a*N+b]);

			fprintf(scalar, "\n");
		}

		fclose(scalar);
	}
	// K:________________________________________________________________________________________
	if(out_scalar==0.3 || out_scalar==0.03 || out_scalar==0.06 || out_scalar==0.006){

		sprintf(name,"../database/K%d.dat", t);
		FILE *scalar = fopen(name,"w");

		for(a=0; a<N; a++){
			for(b=0; b<N; b++)
				fprintf(scalar, "%g ", K[a*N+b]);

			fprintf(scalar, "\n");
		}

		fclose(scalar);
	}
	

}
