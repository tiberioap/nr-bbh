#include "params.h"

void OutVector(int t, double *beta_I, double *GammaTil_I, double *x_I){

	int i, a,b, ab;
	char name[50];

	// beta_I:___________________________________________________________________________
	if(out_vector==1.1 || out_vector==1.32 || out_vector == 1.43 || out_vector == 1.716){

		sprintf(name,"../database/beta_I%d.dat", t);
		FILE *vector = fopen(name,"w");

		for(a=0; a<N; a++){
			for(b=0; b<N; b++){

				ab = a*N+b;
				
				fprintf(vector, "%g\t %g\t %g\t %g\n", a*dx, b*dx,
					beta_I[ab*rank1+0], beta_I[ab*rank1+1]);
			}

			fprintf(vector, "\n");
		}

		fclose(vector);
	}
	// GammaTil_I:_______________________________________________________________________
	if(out_vector==1.2 || out_vector==1.32 || out_vector == 1.56 || out_vector == 1.716){

		sprintf(name,"../database/GammaTil_I%d.dat", t);
		FILE *vector = fopen(name,"w");

		for(a=0; a<N; a++){
			for(b=0; b<N; b++){

				ab = a*N+b;
				
				fprintf(vector, "%g\t %g\t %g\t %g\n", a*dx, b*dx,
					GammaTil_I[ab*rank1+0], GammaTil_I[ab*rank1+1]);
			}

			fprintf(vector, "\n");
		}

		fclose(vector);
	}
	// x_I:______________________________________________________________________________
	if(out_vector==1.3 || out_vector==1.43 || out_vector == 1.56 || out_vector == 1.716){

		FILE *vector = fopen("../database/x_I.dat","a");

		fprintf(vector, "%g\t %g\t %g\t %g\n", x_I[(t*Dm+0)*2+0], 
		x_I[(t*Dm+1)*2+0], x_I[(t*Dm+0)*2+1], x_I[(t*Dm+1)*2+1]);

		fclose(vector);
	}

}
