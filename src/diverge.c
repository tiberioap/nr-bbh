#include "params.h"

void Diverge(int t, int w, double *f, char n[15]){

	if(isinf(f[w])){
 
		printf("\n→  At t = %d, %s is Inf !\n", t, n);
		exit(EXIT_FAILURE);
	}
	
	else if(isnan(f[w])){
 
		printf("\n→  At t = %d, %s is NaN !\n", t, n);
		exit(EXIT_FAILURE);
	}	

	if(fabs(f[w]) > 1e3) \
		printf("At t = %d, %s = %g\n", t, n, f[w]); 

}
