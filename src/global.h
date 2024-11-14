
// Definitions:_______________________________________________________________________________________________________
#define Nt (int)(T/dt)
#define N (int)(L/dx)

#define rank0 (int)(pow(Dm,0))
#define rank1 (int)(pow(Dm,1))
#define rank2 (int)(pow(Dm,2))
#define rank3 (int)(pow(Dm,3))

#define rank0Dm (int)(rank0*pow(N,Dm))
#define rank1Dm (int)(rank1*pow(N,Dm))
#define rank2Dm (int)(rank2*pow(N,Dm))
#define rank3Dm (int)(rank3*pow(N,Dm))

#define r_s (2*M) // Schwarzschild radius (geometrical units)
#define pi acos(-1)

// Prototypes:________________________________________________________________________________________________________
double *Al(int);

void Ic(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

void Diverge(int, int, double *, char n[]);

void OutScalar(int, double *, double *, double *);

void OutVector(int, double *, double *, double *);

void OutTensor(int, double *, double *);

double Del(int, int, int, int, int, double *);

double Del2(int, int, int, int, double *);

double Ddel(int, int, int, int, int, int, double *);

double beta_Del(int, int, int, int, int, double *, double *);

double KO4th(int, int, int, int, double *);

double Lie(int, int, int, int, int, double *, double *);

void Inv(int, int, double *, double *);

void Metric(int, int, int, double *, double *, double *, double *, double *, double *, double *);

void Connections(int, int, double *, double *, double *);

void ConformTensor(int, int, double *, double *, double *, double *);

int epsilon(int, int, int);

double Ricci_chi(int, int, int, int, double *, double *, double *, double *, double *);

double Ricci_til(int, int, int, int, double *, double *, double *, double *, double *);

void CN_scheme(int, int, int, int, int, double *, double *, double *, double *, double *, double *);

double func(int, int, int, int, double *, double *, double *, double *, double *);

