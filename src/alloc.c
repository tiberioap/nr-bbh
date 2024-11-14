#include "params.h"

double *Al(int size){

	double *vector;	

	vector = (double*) malloc (size*sizeof(double));

	if(!vector){
		fprintf(stderr, "Error allocating vector !\n");
		exit(EXIT_FAILURE);
	}

        return vector;
}
