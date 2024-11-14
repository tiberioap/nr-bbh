#include "params.h"

void OutTensor(int t, double *gammaTil_ij, double *ATil_ij){

	int a,b, ab;
	char name[50];

	// gammaTil_ij:_____________________________________________________________________________________
	if(out_tensor==2.1 || out_tensor==4.62){

		sprintf(name,"../database/gammaTil_ij%d.dat", t);
		FILE *tensor = fopen(name,"w");

		for(a=0; a<N; a++){
			for(b=0; b<N; b++)
				fprintf(tensor, "%g\t %g\t %g\n", a*dx, b*dx, gammaTil_ij[(a*N+b)*rank2+0]);

			fprintf(tensor, "\n");
		}

		fclose(tensor);
	}
	// ATil_ij:__________________________________________________________________________________________
	if(out_tensor==2.2 || out_tensor==4.62){

		sprintf(name,"../database/ATil_ij%d.dat", t);
		FILE *tensor = fopen(name,"w");

		for(a=0; a<N; a++){
			for(b=0; b<N; b++)
				fprintf(tensor, "%g\t %g\t %g\n", a*dx, b*dx, ATil_ij[(a*N+b)*rank2+0]);

			fprintf(tensor, "\n");
		}

		fclose(tensor);
	}

}
