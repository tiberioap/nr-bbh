w#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_linalg.h>

//________________________________________________________________
// grid lenght
#define L 200

// time period
#define T 20

#define dx 0.5  
#define dt (0.01*dx) 

// grid dimension
#define Dm 2 

// key to start to compute dynamical eqs
#define key_dynamics 1 // (on →  1, off →  0) 

//__________________________________________________________________
// particles mass
#define M 0.5

#define eta (1./M)

// bh's separation = d*r_s
#define d 15.

#define lapse_Ic 0 // 0 →  = 1.0; 1 →  = (2/psi-1); 2 →  = psi^-2

#define __momentum_spin__ double \
P[] = {0.0,0.0, 0.0,0.0}, \
S[] = {0.0,0.0, 0.0,0.0};

// weigth to compute 1/chi
#define delta 0.3

// Kreiss–Oliger dissipation constant
#define eps 0.001 

//__________________________________________________________________
// number of files
#define Nf 100

// key output: set zero to unable; multiplay the values to combine
#define out_scalar 0.006  // (chi: 0.1, alpha: 0.2, K: 0.3)
#define out_vector 1.2  // (beta_I: 1.1, GammaTil_I: 1.2, x_I: 1.3)
#define out_tensor 0.0  // (gammaTil_ij: 2.1, ATil_ij: 2.2)

#define key_eog 1 // (0 →  off, 1 →  on)

#include "global.h"
