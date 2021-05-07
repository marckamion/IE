void de ( void f ( double t, double y[], double yp[] ), int neqn, double y[], 
  double *t, double tout, double relerr, double abserr, int *iflag, double yy[], 
  double wt[], double p[], double yp[], double ypout[], double phi[], 
  double alpha[], double beta[], double sig[], double v[], double w[], 
  double g[], int *phase1, double psi[], double *x, double *h, double *hold, 
  int *start, double *told, double *delsgn, int *ns, int *nornd, int *k, int *kold, 
  int *isnold );
int i4_sign ( int i );
void intrp ( double x, double y[], double xout, double yout[], double ypout[], 
  int neqn, int kold, double phi[], double psi[] );
void ode ( void f ( double t, double y[], double yp[] ), int neqn, double y[], 
  double *t, double tout, double relerr, double abserr, int *iflag, 
  double work[], int iwork[] );
double r8_abs ( double x );
double r8_add ( double x, double y );
double r8_epsilon ( void );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_sign ( double x );
void step ( double *x, double y[], void f ( double t, double y[], double yp[] ), 
  int neqn, double *h, double *eps, double wt[], int *start, double *hold, 
  int *k, int *kold, int *crash, double phi[], double p[], double yp[], 
  double psi[], double alpha[], double beta[], double sig[], double v[], 
  double w[], double g[], int *phase1, int *ns, int *nornd );
void timestamp ( void );
