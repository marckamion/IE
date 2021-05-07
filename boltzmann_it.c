/*     time evolution of of perturbation variables with Boltzmann eqns  */

/*       now includes neutrino flucts using CLASS UFA   */
/*    TO DO:     
		 figure out allowable/optimal taustep
		 figure out allowable starting time
		 loop over k
                 figure out how to include HyRec?
                 figure out how to use exact evolution of scale factor with tau?  */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ode.h"
#include "nrutil.c"
#include "nrutil.h"
#include "ode.c"
#include "qromb.c"
#include "trapzd.c"
#include "polint.c"
#include "splint.c"
#include "spline.c"
#include "mmid.c"
#include "odeint.c"
#include "rkqc.c"
#include "rk4.c"

#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define PI (3.141592653589793)
#define TINY (0.00000000001)


static double hhh;          /*  hubble parm in units of 100 km/s/Mpc  */
static double Omega_c;      /*   cdm density  */
static double Omega_b;    /*  baryon density  */
static double Omega_r;    /* radiation density today  */
static double Omega_ph;    /*  photon density today   */
static double fnu;   /*   neutrino fraction for N_eff=3.04  */
static double Hinverse;    /*  c/H_0  in Mpc   */
static double tau0,tau_eq,tau_ls;     /* conformal time today, MR-eq, ls   */
static double ne0;       /* electron density today in 1/cm^3  */
static double T0;        /*   CMB temperature today in eV  */
static double rhocrit;   /*  critical density today in GeV/cm^3  */
static double calHeq;    /*  calH ( =a*H )   at equality  */
static double k;         /*   comoving wavenumber in 1/Mpc  */
static double a_eq;      /*    scale factor at MR equality   */
static double taustep,tau_handoff;      /*  integration timestep and handoff from TCA to IE   */
static double *z_array,*xe_array,*z2_array,*xe2_array,*Tm_array,*Tm2_array;   /*  arrays for interpolation of thermo quantities  */
static double sigma;    /*  photon quadrupole  */
static double YHe;         /*   helium fraction   */
static int lmax;



int main(int argc,char *argv[])
{
  double tau0integrand(),xe(),scalefactor(),dotkappa(),calH();
  double j2overx2(),j2(),j2prime(),RLL2(),cs2();
  double tau;   /*  conformal time in Mpc   */
  double kappadot,kappa,visibility;
  double ktauinitial;    /* initial value of k*tau for TCA integration  */
  double tauinitial;
  double tau_end,Pie;
  double kt2,kt3,kt4,om,eps,Rb;
  double *y,*dydt;    /*    TCA perturbation variables  */
  double *workspace;
  int *iworkspace;
  double taumax;    /*  start and end of integration   */
  void odeint(),rkqc(),allderivs();
  void ode(),allderivs();
  int nok,nbad;
  int Npts;     /* number of integration time steps   */
  double *Delta0,*hdot,*alphadot,*thetag,*thetab,*deltab,*deltac,*Delta2,*Pi,*tauj,*kappafrom1,*kappadotprime;
  double R,a,x,RoldcH,thetabprime,dotalpha,kappadd;
  double test1,test2,test3,dekappa,eip1,ei,Deltakappa1,Deltakappa;
  double Deltakappa2,Deltakappa4,Deltakappa3;
  double Kip1,Ki,Kim1,Iip1,Ii,Iim1,Pii,Piip1,Piim1,edk,dk,dk1,thetabpp,sum;
  double x0,x1,x2;
  double z,fraction,Tm;
  double *taucount,*visibilitycount;

  int i,j,count,n;

  FILE *tcafile,*kappadotfile,*testfile,*xefile;

  lmax = 50;

  tcafile=fopen("boltzmannresults.out","w");
  kappadotfile=fopen("kapparesults.out","w");
  testfile=fopen("test.out","w");


  Npts=6000;
  z_array=vector(1,8000);
  xe_array=vector(1,8000);
  z2_array=vector(1,8000);
  xe2_array=vector(1,8000);
  Tm_array=vector(1,8000);
  Tm2_array=vector(1,8000);
  taucount=vector(1,2000);
  visibilitycount=vector(1,2000);


  /*   here are the cosmological parameters; using Planck best-fit values  */
  hhh=0.674;
  Omega_c=0.12038/hhh/hhh;
  Omega_b=0.0224/hhh/hhh;
  fnu=0.6904;
  Omega_ph=2.4735e-5/hhh/hhh;
  Omega_r=Omega_ph*(1.0+fnu);
  Hinverse=2998.0/hhh;
  T0=2.7255 * 8.6217e-5;
  rhocrit = 1.054e-5 * hhh * hhh;      /*  in GeV/cm^3   */
  ne0 = rhocrit*Omega_b/0.938;    /*   this is actually the *nucleon* density  */
  YHe = 0.245;

  a_eq = Omega_r/(Omega_c + Omega_b);     /*  scale factor at MR equality  */

  /* calculate conformal time today, at MR equality, and at ls  */
  tau0=qromb(tau0integrand,0.0000000000001,1.0) * Hinverse;
  tau_eq=qromb(tau0integrand,0.0000000000001,Omega_r/(Omega_c+Omega_b)) * Hinverse;
  tau_ls=qromb(tau0integrand,0.0000000000001,1.0/1100.0) * Hinverse;
  /*   diagnostic  print  */
  /*printf("tau0=%f     tau_eq=%f  tau_ls=%f  z_eq=%f \n",tau0,tau_eq,tau_ls,1.0/(Omega_r/(Omega_c+Omega_b)));  */
  calHeq=calH(tau_eq);     /* calH is Hubble paramter times scale factor  */

  /*  set up lookup table to evaluate xe(tau)  from output of HyRec-2  */
  xefile = fopen("output_xe.dat","r");
  n=8000;
  while(fscanf(xefile,"%lf %lf %lf",&z,&fraction,&Tm) != EOF) {
    z_array[n]=z;
    xe_array[n]=fraction;
    Tm_array[n]=Tm;
    n--;
  }
  fclose(xefile);
  spline(z_array,xe_array,8000,1.0e40,1.0e40,xe2_array);
  spline(z_array,Tm_array,8000,1.0e40,1.0e40,Tm2_array);
  /*   print out kappa, dotkappa, xe  for check   */
  kappa=0.0;   count=0;
  for(tau = tau_ls*5.0; tau>=tau_ls/20.0; tau -= (5.0-1.0/20.0)*tau_ls/1000.) {
    count++;
    taucount[count]=tau;
    kappadot=dotkappa(tau);
    kappa += kappadot * 4.5*tau_ls/1000.0;
    visibilitycount[count] = kappadot*exp(-kappa);
  }
  while(count>0) {
    fprintf(kappadotfile,"%f   %f   %f   %f   %f  %e\n",taucount[count],1.0/scalefactor(taucount[count]),xe(taucount[count]), dotkappa(taucount[count]),visibilitycount[count],cs2(taucount[count]));
    count = count-1;
  }
  fclose(kappadotfile);

  y=vector(1,9 + 2*(lmax-1));      /*   number of perturbation variables  */
  dydt=vector(1,9+2*(lmax-1));      /*   number of perturbation variables  */

  workspace=vector(1,101+21*(9+2*(lmax-1)));
  iworkspace=ivector(1,5);


  /*   loop over k values  */
  /*  for(k=5.0/tau0; k<= 3001.0/tau0; k+=30.0/tau0) {  */
  for(k=0.2; k<0.22; k+=1.0) {

    tauinitial=10.0;       /*  start of TCA  integration  */
    tau_handoff=160.0;     /* handoff from TCA to IE  */
    tau_end = 450.0;       /*   end of integration  */
    taustep = 0.1;      /*  size of conformal-time steps  */


    /*   y[1] = hdot   y[2] = Delta0    y[3] =  thetag   y[4] = deltac  y[5] = deltaab  y[6]=thetab     y[7] = delta_nu   y[8]=theta_nu   y[9]=sigma_nu   */
    /*   y[10]=Delta_2  y[11]=Delta_3 ...... y[9+(lmax-1)] = Delta_lmax 
	 y[9+lmax] = DeltaE_2 ..... y[9+2(lmax-1)] =DeltaE_lmax  */

    /*  to simplify prefactors, the DeltaE_l here will be those without the 
	(3/4) sqrt((l+2)!/(l-2)!)  and Pi will be DeltaT2 + 9*DeltaE2   */

    /*   initial conditions   from Cyr-Racine & Sigurdson; this is probably overkill but might be enable TCA to be started a bit later, especially for smaller k  */
    ktauinitial = k*tauinitial;
    kt2=ktauinitial*ktauinitial;
    kt3=ktauinitial*kt2;
    kt4 = kt2*kt2;
    om = (1.0/Hinverse)*(Omega_b+Omega_c)/sqrt(Omega_r);
    eps = 1.0/dotkappa(tauinitial)/tauinitial;
    Rb=Omega_b/(Omega_c+Omega_b);
  
    y[1] = 2.0 * k * ktauinitial - (3.0/5.0)*om*kt2 - (1.0/9.0)*(-10.0/(4.0*fnu+15.0) + 1.0)*kt3*k + (1.0/4.0)*om*om*kt2*tauinitial;
    y[2] = -kt2/6.0 * (1.0-om*tauinitial/5.0) +(4.0*fnu+10.0)/4.0/27.0/(4.0*fnu+15.0)*kt4 - om*om*kt2*tauinitial*tauinitial/96.0;
    y[3] = -( k/18.0 * kt3 + 8.0*k*kt3*eps/(36.0*fnu+135.0) - om*kt4*(1.0+5.0*Rb-fnu)/120.0/(1.0-fnu) + 2.0*(2.0*(5.0*Rb-9.0)*fnu+75.0*Rb+8.0*fnu*fnu+10.0)/15.0/(fnu-1.0)/(2.0*fnu+15.0)/(4.0*fnu+15.0)*om*kt4*eps + 16.0*kt3*k*eps*(6.0*fnu+181.0)/45.0/(2.0*fnu+15.0)/(4.0*fnu+15.0));
    y[4] = - kt2/2.0 +0.1*om*kt2*tauinitial+kt4/72.0*(-10.0/(4*fnu+15.0)+1.0) - om*om*kt2*tauinitial*tauinitial/32.0;
    y[5] = (3.0/4.0)*y[2];
    y[6] = y[3];
    y[7] = -(2.0/3.0)*kt2 * (1.0-om*tauinitial/5.0);
    y[8] = -(1.0/18.0)*(4.0*fnu+23.0)/(4.0*fnu+15.0)*kt3*k;
    y[9] = 4.0/(12.0*fnu+45.0)*kt2;
    for(n=10;n<=9+2*(lmax-1);n++) {
      y[n]=0.0;
      dydt[n]=0.0;
    }
      
  
    kappadot=dotkappa(tauinitial);

    count=0;      /*  this counts the timestep   */
    kappa=0.0;    /*  we'll measure kappa from the start of the integration  */
    for(tau=tauinitial; tau<=tau_end; tau+=taustep) {

      /* odeint(y,9+2*(lmax-1),tau,tau+taustep,1.0e-3,taustep/10.0,0.0,&nok,&nbad,allderivs,rkqc);  */
	nok=1;
	ode(allderivs,9+2*(lmax-1),y,&tau,tau+taustep,1.0e-5,1.0e-4,&nok,workspace,iworkspace);
	tau=tau-taustep;

      Pie = y[10]+9.0*y[9+lmax];
      a=scalefactor(tau+taustep)/a_eq;
      R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* the inverse of MB' and CLASS's  R   */
      dotalpha = y[1] - 3.0 * sqr(calHeq)/a/a* (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);

      fprintf(tcafile,"%d  %f  %e  %e  %e  %e  %e  %e\n",0,tau+taustep,y[1],y[2],y[6],dotalpha,2.0*y[10],Pie);
    }
  }

	      /*     fprintf(tcafile,"%1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.6e  %1.6e  %1.3e  %1.3e\n",k,tau+taustep,y[1],y[2],y[3],y[4],y[5],y[6],2.0*y[10],Pie,y[3]-y[6],
		     (1.0/k/k/k)*(Omega_c *y[4]+Omega_b*y[5])*(Omega_c *y[4]+Omega_b*y[5]) );  */


  fclose(tcafile);
  fclose(testfile);
}
  

double xe(double tau)    /*   ionization fraction from HyRec-2 lookup table  */
{                        /*  this is actually n_e/n_p    */
  void splint();
  double scalefactor();
  double y,redshift;

  redshift = 1.0/scalefactor(tau) - 1.0;

  if(redshift<8000.) {
    splint(z_array,xe_array,xe2_array,8000,redshift,&y);
  }
  else {
    y = xe_array[8000];
  }

  return y;
}


double Tb(double tau)    /*   ionization fraction from HyRec-2 lookup table  */
{                        /*  this is actually n_e/n_p    */
  void splint();
  double scalefactor();
  double y,redshift;

  redshift = 1.0/scalefactor(tau) - 1.0;

  if(redshift<8000.0) {
    splint(z_array,Tm_array,Tm2_array,8000,redshift,&y);
  }
  else {
    y = Tm_array[8000];
  }

  return y;
}


double cs2(double tau)     /*   return sound speed squared   */
{
  void splint();
  double scalefactor(),Tb(),xe();
  double y,redshift(),mu,dlnTdlna;

  mu = 0.938 / ( (1.0-YHe)*(1.0+xe(tau))+YHe/4.0);

  dlnTdlna = scalefactor(tau)/Tb(tau) * ( Tb(1.01*tau)-Tb(0.99*tau)) / (scalefactor(1.01*tau)-scalefactor(0.99*tau));

  return Tb(tau) * 8.617e-14/mu * (1.0 - (1.0/3.0)*dlnTdlna);

  /*     8.617e-14 GeV/K     mu = mp * (np + 4 * nHe )/( np*(1+xe)+ nHe)  */
}




double scalefactor(double tau)     /*  use early-time approximation;  returns a/e_eq */
{
  double sqrt2m1=0.4142135623731;
  double y;

  y=  sqrt2m1 * (tau/tau_eq);

  return Omega_r/(Omega_c+Omega_b) * y*(y+2.0);
}


double calH(double tau)    /*  \dot a/a   */
{
  double sqrt2m1=0.4142135623731;
  double y;

  y=  sqrt2m1 * (tau/tau_eq);

  return 2.0 * (1.0+y) / tau / (2.0+y);
}



/*   i don't think this is used anywhere  */  
double conformaltime(double a)
{
  double sqrt2m1=0.4142135623731;

  return tau_eq * (sqrt( (a/a_eq)+1.0) - 1.0)/sqrt2m1;
}
  

double dotkappa(double tau)    /*   \dot\kappa in units of 1/Mpc  */
{
  double xe(),a,scalefactor();

  a=scalefactor(tau);

  return 2.052 * ne0 * (1.0-YHe) * xe(tau) / a/a;    /* 2.052 = thomson cross section in cm^2  times    Mpc/cm  */
}
  


/*  integrand to determine conformal time as a function of scale factor  */  
double tau0integrand(double a)
{
  double E();

  return 1.0/a/a/E(1.0/a-1.0);
}

  
double E(double z)
{
  return sqrt( (Omega_b+Omega_c)*pow(1.0+z,3.0)+Omega_r*pow(1.0+z,4.0)+(1.0-Omega_b-Omega_c-Omega_r));
}


void allderivs(double t, double *y, double *dydt)
/*    differential equations for TCA and full Boltzmann hierarchy   */
{
  double calH(),a,R,cH,dotkappa(),kappadot,slip,dotalpha,kappadd,Pi,ell,thetabprime;
  double aprime_over_a,cHprime,thetadot0,alphadd,F,Fdot,sigmadot,slip1,sigma1;
  double ss2,cs2();
  int l;
  
  /*   y[1] = hdot   y[2] = Delta0    y[3] =  thetag   y[4] = deltac  y[5] = deltaab  y[6]=thetab   */

  ss2=cs2(t);     /*    get the sound speed squared  */
  a=scalefactor(t)/a_eq;
  R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* note that this is the inverse of MB' and CLASS's  R   */
  cH = calH(t);
  cHprime = (calH(t*1.01)-calH(0.99* t))/(0.02*t);
  kappadot=dotkappa(t);

  if(t>=tau_handoff) {
    dotalpha = y[1] - 3.0 * sqr(calHeq)/a/a * (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);
  
    dydt[1]= - cH*y[1] -3.0 * sqr(calHeq) * ( (4.0/a/a) * (Omega_ph/Omega_r) * ( y[2] + fnu*y[7]/4.0) + 0.5 * Omega_c/(Omega_b+Omega_c) /a * y[4] +0.5*y[5]*Omega_b/(Omega_b+Omega_c)/a);
    dydt[2] = -y[3]/3.0 - y[1]/6.0;
    dydt[6] = -cH*y[6]+ ss2*k*k*y[5] + kappadot/R*(y[3]-y[6]);
    dydt[3] =  -R*(dydt[6]+cH*y[6] - ss2*k*k*y[5])+ k*k*(y[2]-2.0*y[10]);
    dydt[4] = -y[1]/2.0;
    dydt[5] = -y[6]-y[1]/2.0;
    dydt[7] = -(4.0/3.0)*y[8] - (2.0/3.0)*y[1];
    dydt[8] = k*k*(y[7]/4.0 - y[9]);
    dydt[9] = -(3.0/t) * y[9] + (2.0/3.0)*y[8] + (1.0/3.0)* y[1];

    Pi = y[10]+9.0*y[9+lmax];
    dydt[10] = dotalpha/15.0 - kappadot*(y[10]- Pi/10.0) + (2.0/5.0)*y[3]/3.0 - (3.0*k/5.0)*y[11];  /*   photon quadrupole  */
    dydt[9+lmax] = -kappadot*y[9+lmax]+(1.0/15.0)*kappadot*Pi - k*y[9+lmax+1]; /*  polarization quadrupole  */
    for(l=3;l<=lmax-1;l++) {
      ell = (double)l;
      dydt[l+8] = -kappadot*y[l+8] + ell*k/(2.0*ell+1.0)*y[l+7] - (ell+1.0)*k/(2.0*ell+1.0)*y[l+9];
      dydt[l+7+lmax] = -kappadot*y[l+7+lmax] + (ell-2.0)*k/(2.0*ell+1.0) * y[l+6+lmax] - (ell+3.0)*k/(2.0*ell+1.0) * y[l+8+lmax];
    }
    ell=(double)lmax;
    dydt[lmax+8] = k*y[lmax+7] - ((ell+1.0)/t + kappadot)*y[lmax+8];
    dydt[2*lmax+8] = k*y[2*lmax+7] - ((ell+1.0)/t + kappadot)*y[2*lmax+8];
  }

  else {
     dydt[1]= - cH*y[1] -3.0 * sqr(calHeq) * ( (4.0/a/a) * (Omega_ph/Omega_r) * ( y[2] + fnu*y[7]/4.0) + 0.5 * Omega_c/(Omega_b+Omega_c) /a * y[4] +0.5*y[5]*Omega_b/(Omega_b+Omega_c)/a);
  dydt[8] = k*k*(y[7]/4.0 - y[9]);

  /*   use TCA for sigma at early times  */
    dotalpha = y[1] - 3.0 * sqr(calHeq)/a/a * (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);
    
    aprime_over_a = cH*cH+cHprime;
    kappadd = (dotkappa(1.01*t)-dotkappa(0.99*t))/(0.02*t);
    thetadot0 = R/(1.0+R) *(-cH*y[3] + ss2*k*k*y[5]+k*k/R*y[2]);

    alphadd = dydt[1] - 2.0*cH * (dotalpha-y[1]) - 3.0*sqr(calHeq)/a/a * Omega_ph/Omega_r * (thetadot0*(1.0+R) + cH*R*y[3] + fnu * dydt[8]);
    
  /*  lowest-order TCA approx for sigma (photon quadrupole) from CLASS paper  */
    sigma = 8.0/45.0/kappadot*(2.0 * y[3]+ dotalpha);
    sigma1= sigma;

    /*  now get sigma to second order in 1/kappadot  */
    sigma = (1.0+11.0/6.0*kappadd/kappadot/kappadot)*sigma - (11.0/6.0)/kappadot *  (8.0/45.0)/kappadot* (2.0*thetadot0+alphadd);

    F = R/(1.0+R)/kappadot;

    slip = -(kappadd/kappadot + 2.0*cH*R/(1.0+R))*(y[3]-y[6]) - F * (-aprime_over_a * y[6] + k*k*(-2.0*cH*y[2] + ss2*(-y[6]-y[1]/2.0)- ( -y[3]/3.0-y[1]/6.0)));
    slip1 = slip;

    /*  now get slip to partial 2nd order (compromise_class approx)  */
    sigmadot = -kappadd/kappadot*sigma+(8.0/45.0)/kappadot*(2.0*thetadot0+alphadd);
    Fdot= 1.0/(1.0+R)/(1.0+R)*cH*R -R/(1.0+R)*kappadd/kappadot/kappadot;
    slip = (1.0- 2.0*cH*F)*slip - F*k*k*(2.0*cH*sigma+sigmadot-(1.0/3.0 - ss2)*(F*thetadot0 +  2.0* Fdot*y[6]));

    /*  and here's the TCA approx for thetabdot  */
    thetabprime= - R /(1.0+R) * (cH*y[6]-ss2*k*k*y[5] - k*k/R*(y[2]-sigma)+slip/R);
  dydt[2] = -y[3]/3.0 - y[1]/6.0;
  dydt[3] =  -R*(thetabprime+cH*y[6] - ss2*k*k*y[5])+ k*k*(y[2]-sigma);
  dydt[4] = -y[1]/2.0;
  dydt[5] = -y[6]-y[1]/2.0;
  dydt[6] = thetabprime;
  dydt[7] = -(4.0/3.0)*y[8] - (2.0/3.0)*y[1];
  dydt[9] = -(3.0/t) * y[9] + (2.0/3.0)*y[8] + (1.0/3.0)* y[1];


  for(l=10;l<=9+2*(lmax-1); l++ ) {
    dydt[l]=0.0;
  }
  y[10]=sigma/2.0;
  y[9+lmax]= ((3./2.)*y[10] - (1.0/3.0)*alphadd/kappadot - (2.0/3.0)*thetadot0/kappadot)/9.0;
  }
}






double j2overx2(double x) {
  double x2,x4;

  x2=x*x;
  x4=x2*x2;
  
  if(x>=0.01) {
    return ((-3.0 * x * cos(x) - (x2-3.0)*sin(x))/x2/x)/x2;
  }
  else {
    return 1.0/15.0 - x2/210.0 + x4/7560.0;
  }
}


double j2(double x) {
  double x2,x4;

  x2=x*x;
  x4=x2*x2;
  
  if(x>=0.01) {
    return ((-3.0 * x * cos(x) - (x2-3.0)*sin(x))/x2/x);
  }
  else {
    return 1.0/15.0*x2 - x4/210.0 + x2*x4/7560.0;
  }
}


double j2prime(double x) {
  double x2,x4;

  x2=x*x;
  x4=x2*x2;

  if(x>0.01) {
    return (-x * (x2-9.0) * cos(x) + (4 * x2 -9.0) * sin(x))/x4;
  }
  else {
    return x * (2.0/15.0 - 2.0*x2/105.0 + x4/1260.0);
  }
}


double RLL2(double x) {
  double x2,x4;

  x2=x*x;
  x4=x2*x2;

  if(x>0.1) {
    return - ( 6.0*x*(x2-9.0)*cos(x) + (54.0-24.0*x2+x4)*sin(x))/x/x4;
  }
  else {
    return -1.0/5.0 + 11.0*x2/210.0 - x4/280.0 + 17.0/166320.0*x4*x2;
  }
}


  



/*
#include "stiff.c"
#include "stifbs.c"
#include "jacobn.c"
#include "/Users/kamion/SCIENCE/recipes_c/lubksb.c"
#include "/Users/kamion/SCIENCE/recipes_c/ludcmp.c"
#include "simpr.c"
#include "/Users/kamion/SCIENCE/recipes_c/pzextr.c"
#include "/Users/kamion/SCIENCE/recipes_c/odeint.c"
#include "/Users/kamion/SCIENCE/recipes_c/bsstep.c"
#include "rkqc.c"
#include "/Users/kamion/SCIENCE/recipes_c/rk4.c"
#include "/Users/kamion/SCIENCE/recipes_c/rzextr.c"
*/



/*
  free_vector(y,1,9+2*(lmax+1));
  free_vector(dydt,1,9+2*(lmax+1));


  free_vector(z_array,1,8000);
  free_vector(xe_array,1,8000);
  free_vector(z2_array,1,8000);
  free_vector(xe2_array,1,8000);
  free_vector(Tm_array,1,8000);
  free_vector(Tm2_array,1,8000);
  free_vector(taucount,1,2000);
  free_vector(visibilitycount,1,2000);
*/
