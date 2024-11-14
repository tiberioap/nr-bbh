#include "params.h"

double Del(int rank, int l, int ij, int a, int b, double *f){

	// l   →  del index
	// ij  →  tensor index: if rank = 0, ij = 0
	// a,b →  grid indexes

	int n = pow(Dm,rank); // number of coefficients 	
	int ap,am, bp, bm; // (p)lus and (m)inus indexes
	int wp, wm, wpp, wmm;

	double del;

	if(l == 0){

		ap = (a+1)%N;	am = (a-1+N)%N;
		wp = ap*N+b;	wm = am*N+b;
	}
	else if(l == 1){

		bp = (b+1)%N;	bm = (b-1+N)%N;
		wp = a*N+bp;	wm = a*N+bm;
	}
	else
		printf("\nError in del.c: Del() !\n");


	del = (0.5/dx)*(f[wp*n+ij]-f[wm*n+ij]);
	
	return del;	
}

double Del2(int rank, int ij, int a, int b, double *f){

	// if rank = 0, then ij = 0

	int n = pow(Dm,rank); 	
	int ap,am, bp, bm; 
	int ab, apb,amb, abp,abm;
	
	double d1, d2;

	ab = a*N+b;

	ap = (a+1)%N;   am = (a-1+N)%N;
	bp = (b+1)%N;   bm = (b-1+N)%N;

	apb = ap*N+b;   amb = am*N+b;
	abp = a*N+bp;   abm = a*N+bm;

	d1 = pow(dx,-2)*(f[apb*n+ij]-2*f[ab*n+ij]+f[amb*n+ij]);
	
	d2 = pow(dx,-2)*(f[abp*n+ij]-2*f[ab*n+ij]+f[abm*n+ij]);

	if(Dm == 3) // 3-dimensional grid
		printf("\nError in del.c: Del2() !\n");


	return d1 + d2;
}

double Ddel(int rank, int i, int j, int kl, int a, int b, double *f){

	// if rank = 0, then kl = 0

	int n = pow(Dm,rank);  	
	int ap,am, bp,bm; 
	int ab, apb,amb, abp,abm;
	int apbp,ambp, apbm,ambm;
	double ddel;

	ab = a*N+b;

	if(i == j){

		if(i == 0){

			ap = (a+1)%N;	am = (a-1+N)%N;
			apb = ap*N+b;	amb = am*N+b;
	
			ddel = pow(dx,-2)*(f[apb*n+kl]-2*f[ab*n+kl]+f[amb*n+kl]);
		}
		else if(i == 1){
			
			bp = (b+1)%N;	bm = (b-1+N)%N;
			abp = a*N+bp;	abm = a*N+bm;
	
			ddel = pow(dx,-2)*(f[abp*n+kl]-2*f[ab*n+kl]+f[abm*n+kl]);
		}
	}

	else if((i==0 && j==1) || (i==1 && j==0)){
		
		ap = (a+1)%N;	am = (a-1+N)%N;
		bp = (b+1)%N;	bm = (b-1+N)%N;

		apbp = ap*N+bp;  ambp = am*N+bp;	
		apbm = ap*N+bm;  ambm = am*N+bm;	

		ddel = 0.25*pow(dx,-2)*(f[apbp*n+kl]-f[apbm*n+kl]-f[ambp*n+kl]+f[ambm*n+kl]); 
	}
	else
		printf("\nError in del.c: Ddel() !\n");
	

	return ddel;
}

double beta_Del(int rank, int l, int ij, int a, int b, double *beta_I, double *f){

	// if rank = 0, then ij = 0

	int ab = a*N+b;
	int n = pow(Dm,rank);  	
	int ap1,ap2,ap3,am1,am2,am3;
	int bp1,bp2,bp3,bm1,bm2,bm3; 
	int wp1,wp2,wp3, wm1,wm2,wm3;

	double del, beta_del;

	if(beta_I[ab*rank1+l] != 0.0){

		ap1 = (a+1)%N;	am1 = (a-1+N)%N;
		bp1 = (b+1)%N;	bm1 = (b-1+N)%N;

		ap2 = (a+2)%N;	am2 = (a-2+N)%N;
		bp2 = (b+2)%N;	bm2 = (b-2+N)%N;
		
		ap3 = (a+3)%N;	am3 = (a-3+N)%N;
		bp3 = (b+3)%N;	bm3 = (b-3+N)%N;
		
		if(l == 0){
		
			wp1 = ap1*N+b;  wp2 = ap2*N+b;  wp3 = ap3*N+b;
			wm1 = am1*N+b;  wm2 = am2*N+b;  wm3 = am3*N+b;
		}		

		else if(l == 1){
		
			wp1 = a*N+bp1;  wp2 = a*N+bp2;  wp3 = a*N+bp3;
			wm1 = a*N+bm1;  wm2 = a*N+bm2;  wm3 = a*N+bm3;
		}		

		else 
			printf("\nError in del.c: beta_Del() !\n");


		if(beta_I[ab*rank1+l] > 0)
			del = f[wp3*n+ij] - 6*f[wp2*n+ij] + 18*f[wp1*n+ij] - 10*f[ab*n+ij] - 3*f[wm1*n+ij];	
				
		else if(beta_I[ab*rank1+l] < 0)
			del = -f[wm3*n+ij] + 6*f[wm2*n+ij] - 18*f[wm1*n+ij] + 10*f[ab*n+ij] + 3*f[wp1*n+ij];	
				

		beta_del = (1/12.)*(beta_I[ab*rank1+l]/dx)*del;
	}

	else beta_del = 0.0;	

	return beta_del;
}

double KO4th(int rank, int ij, int a, int b, double *f){

	// if rank = 0, then ij = 0

	int n = pow(Dm,rank);
	int ab = a*N+b;
	int ap1,am1, bp1, bm1,
	    ap2,am2, bp2, bm2;
	int ap1b,am1b,abp1,abm1,
	    ap2b,am2b,abp2,abm2;

	double d1, d2;

	ap1 = (a+1)%N;   am1 = (a-1+N)%N;
	ap2 = (a+2)%N;   am2 = (a-2+N)%N;

	bp1 = (b+1)%N;   bm1 = (b-1+N)%N;
	bp2 = (b+2)%N;   bm2 = (b-2+N)%N;

	ap1b = ap1*N+b;  am1b = am1*N+b;  abp1 = a*N+bp1;  abm1 = a*N+bm1;
	ap2b = ap2*N+b;  am2b = am2*N+b;  abp2 = a*N+bp2;  abm2 = a*N+bm2;

	d1 = f[ap2b*n+ij]-4*f[ap1b*n+ij]+6*f[ab*n+ij]-4*f[am1b*n+ij]+f[am2b*n+ij];      
	d2 = f[abp2*n+ij]-4*f[abp1*n+ij]+6*f[ab*n+ij]-4*f[abm1*n+ij]+f[abm2*n+ij];      

	return -(eps/dx)*(d1 + d2);
}
