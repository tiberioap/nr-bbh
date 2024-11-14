#include "params.h"

void Inv(int t, int ab, double *A_ij, double *A_IJ){

	int i,j, ij, s;
	double C[rank2], invC[rank2]; // auxiliary matrix

	gsl_matrix_view m;
	gsl_matrix_view invm;

	gsl_permutation *p = gsl_permutation_alloc(Dm);

	//__________________________________________________
	for(i=0; i<Dm; i++)
		for(j=0; j<Dm; j++){
	
			ij = i*Dm+j;

			C[ij] = A_ij[ab*rank2+ij];
		}

	//__________________________________________________
	m = gsl_matrix_view_array(C,Dm,Dm);

	invm = gsl_matrix_view_array(invC,Dm,Dm);

	// computing inverse matrix ________________________
	gsl_linalg_LU_decomp(&m.matrix, p, &s);    

	gsl_linalg_LU_invert(&m.matrix, p, &invm.matrix);

	//__________________________________________________
	for(i=0; i<Dm; i++)
		for(j=0; j<Dm; j++){
	
			ij = i*Dm+j;
	
			A_IJ[ab*rank2+ij] = invC[ij];
		}


	gsl_permutation_free(p);
}
