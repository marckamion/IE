/*     This code calculates the time evolution of cosmological-perturbation 
       variables by solving integral equations for the photon-intensity
       quadrupole (Delta2) and photon-polarization quadrupole (Pi).  The
       other quantities are evolved with differential equations as in the
       standard approach.
       This version considers only one k value (k=0.2) that correspond roughly to 
       a CMB  l ~ 3000.  Neutrinos are assumed to
       be massless and are treated as a GDM with w=c_s^2=c_vis^2=1/3
       The integration stops at a conformal time tau=450 Mpc, after 
       recombination but well before reionization.

       Perturbation variables are output as a function of conformal time 

       !!!!THIS CODE IS NOT CERTIFIED FOR SCIENCE!!!.....although the
       warning is not required given that the code stops short of
       calculating CMB power spectra

       The ionization history and baryon temperature are read in from
       the output of HyRec-2, here assumed to be in the file output_xe.dat

       The main results are output to the file iterative.out.  The quantities
       exported are  tau, hdot, Delta0, alphadot, thetab, Delta2, and Pi 
       (i.e., those that appear as sources in the CMB transfer functions

       Other things:
       - this code uses an analytic approximation for the early-time expansion
       history.  

       - This code uses Planck best-fit cosmological parameters
*/


/*   standard include files, and some for integration, interpolation, and ODEs  */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.c"
#include "qromb.c"
#include "ode.c"
#include "ode.h"
#include "trapzd.c"
#include "polint.c"
#include "splint.c"
#include "spline.c"
#include "odeint.c"
#include "rkqc.c"
#include "rk4.c"
#include "mmid.c"

#define sqr(x) ((x)*(x))

/************************************************************************************/

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
static double tau_handoff;      /*  integration timestep and handoff from TCA to IE   */
static double *z_array,*xe_array,*z2_array,*xe2_array,*Tm_array,*Tm2_array;   /*  arrays for interpolation of thermo quantities  */
static double *taugrid,*Delta2,*Delta22_array,*Pi,*Pi2_array;
static double sigma,sigmaold,Pie;    /*  photon quadrupoles  */
static double YHe;         /*   helium fraction   */
static int Npts;     /* number of integration time steps   */

int main(int argc,char *argv[])
{
  double tau0integrand(),xe(),scalefactor(),dotkappa(),calH(),rsintegrand();
  double j2overx2(),j2(),j2prime(),RLL2(),cs2();
  double tau,DeltaTau;   /*  conformal time in Mpc  and time step for TCA */
  double kappadot,kappa,visibility;    
  double ktauinitial;    /* initial value of k*tau for TCA integration  */
  double tauinitial;     /* initial tau for the TCA   */
  double tau_end;    /*  final tau value  */
  double kt2,kt3,kt4,om,eps,Rb;    /*   for initial conditions  */
  double *y,*dydt,*ylate;    /*    all perturbation variables except Delta2 and Pi  */
  void ode(),odeint(),rkqc(),tcaderivs(),exactderivs();    /*  integrators   */
  int nok,nbad;    /*  integer inputs for ODE solvers  */
  double *Delta0,*hdot,*alphadot,*thetag,*thetab,*deltab,*deltac,*kappafrom1,*kappadotprime;   /* arrays to store variables  */
  double R,a,x,Rold,cH,thetabprime,dotalpha,kappadd;
  double test1,test2,test3,dekappa,eip1,ei,Deltakappa1,Deltakappa;
  double Deltakappa2,Deltakappa4,Deltakappa3;
  double Kip1,Ki,Kim1,Iip1,Ii,Iim1,Pii,Piip1,Piim1,edk,dk,dk1,thetabpp,sum;
  double x0,x1,x2;
  double z,fraction,Tm;
  double *taucount,*visibilitycount;
  double *workspace,taulow,tauhigh;
  int *iworkspace;
  double *r_s;
  double rsmatch;
  double Delta2_interp(),Pi_interp();
  double aa,bb,cc,dd,D2,ppp;

  int i,j,count,n,imax,iteration;
  int swit,lowercount,startcount;
  swit=0;

  FILE *tcafile,*kappadotfile,*xefile,*testfile;
  y=vector(1,9);      /*   number of TCA perturbation variables  */
  ylate=vector(1,9);      /*   number of TCA perturbation variables  */
  dydt=vector(1,9);      /*   number of TCA perturbation variables  */

  workspace=vector(1,101+21*9);
  iworkspace=ivector(1,5);

  Npts=401;      /*   (less than) the number of tau grid points  */
  /****   These are arrays of the perturbation variables    */
  taugrid=vector(1,Npts);      /*   grid of tau points  */
  r_s=vector(1,Npts);       /*   sound horizon    */
  Delta0=vector(1,Npts);
  Delta2=vector(1,Npts);
  Delta22_array=vector(1,Npts);
  Pi2_array=vector(1,Npts);
  hdot=vector(1,Npts);
  alphadot=vector(1,Npts);
  deltac=vector(1,Npts);
  deltab=vector(1,Npts);
  thetag=vector(1,Npts);
  thetab=vector(1,Npts);
  Pi=vector(1,Npts);
  kappafrom1=vector(1,Npts);
  kappadotprime=vector(1,Npts);

  /*  These are arrays to read in output of HyRec-2   */
  z_array=vector(1,8000);
  xe_array=vector(1,8000);
  z2_array=vector(1,8000);
  xe2_array=vector(1,8000);
  Tm_array=vector(1,8000);
  Tm2_array=vector(1,8000);

  /*  these are arrays to evaluate arrays for the visibility function  */
  taucount=vector(1,2000);
  visibilitycount=vector(1,2000);

  /*   here are the cosmological parameters; using Planck best-fit values  */
  hhh=0.674;
  Omega_c=0.12038/hhh/hhh;
  Omega_b=0.0224/hhh/hhh;
  /* Omega_b = 0.021952/hhh/hhh;  */
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
  /*   print out kappa, dotkappa, xe  for check

       kappadotfile=fopen("kapparesults.out","w");
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
  */


  /*  initialize the tau and the grid and the values of kappa and r_s on it  */
  i=1;
  tau=160.0;
  taugrid[i]=tau;
  kappafrom1[i]=0.0;
  kappadotprime[i]=dotkappa(tau);
  r_s[i]=qromb(rsintegrand,0.0,tau);
  while(tau<=240.0) {
    i++;
    tau+=1.0;    /*   was 2.5  for Npts=182  and  1.0 for Npts=401   */
    taugrid[i]=tau;
    kappafrom1[i]=kappafrom1[i-1] +  qromb(dotkappa,taugrid[i-1],tau);
    kappadotprime[i]=dotkappa(tau);
    r_s[i] = r_s[i-1]+qromb(rsintegrand,taugrid[i-1],tau);
  }
  while(tau<=350.0) {
    i++;
    tau+=0.5;      /*    was 1.0 for Npts=182 and 0.5 for Npts=401 */
    taugrid[i]=tau;
    kappafrom1[i]=kappafrom1[i-1] +  qromb(dotkappa,taugrid[i-1],tau);
    kappadotprime[i]=dotkappa(tau);
    r_s[i] = r_s[i-1]+qromb(rsintegrand,taugrid[i-1],tau);
  }
  while(tau<=450.0) {
    i++;
    tau+=1.0;       /*    was 2.5 for Npts=182 and 1.0 for Npts=401 */
    taugrid[i]=tau;
    kappafrom1[i]=kappafrom1[i-1] +  qromb(dotkappa,taugrid[i-1],tau);
    kappadotprime[i]=dotkappa(tau);
    r_s[i] = r_s[i-1]+qromb(rsintegrand,taugrid[i-1],tau);
  }
  imax=i;
  if(imax != Npts) {
    printf("imax != Npts;   imax=%d\n",imax);
    exit(1);
  }
      
  /*
    testfile=fopen("test.out","w");
    for(i=1;i<=imax;i++) {
    fprintf(testfile,"%d  %f  %f  %f\n",i,taugrid[i],kappafrom1[i],r_s[i]);
    }
    fclose(testfile);
    exit(1);
  */


  k=0.2;       /*   set k value  */

  /*****************************************************************************************/
  /* do tight coupling approximation until tau=160.0   */
  
  tauinitial=10.0;       /*  start of TCA  integration  */
  tau_handoff=160.0;     /* handoff from TCA to IE  */
  tau_end = 450.0;       /*   end of integration  */


  /*  Solve TCA equations from tau = 10.0 to tau=160.0   */
  /*   y[1] = hdot   y[2] = Delta0    y[3] =  thetag   y[4] = deltac  y[5] = deltaab  y[6]=thetab     y[7] = delta_nu   y[8]=theta_nu   y[9]=sigma_nu   */

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
  
  kappadot=dotkappa(tauinitial);

  count=0;      /*  this counts the timestep   */
  kappa=0.0;    /*  we'll measure kappa from the start of the integration  */
  sigmaold=0.0;

  /*   integrate tight-coupling eqns from tauinitial to tau_handoff; do so via steps
        so we can make plots   */
  tcafile=fopen("iterative.out","w");

    DeltaTau = -(tauinitial-tau_handoff)/100.0;
    for(tau=tauinitial; tau<=tau_handoff-0.9*DeltaTau; tau+=DeltaTau) {
      odeint(y,9,tau,tau+DeltaTau,1.0e-6,DeltaTau/10.0,0.0,&nok,&nbad,tcaderivs,rkqc);

      tcaderivs(tau+DeltaTau,y,dydt);
      a=scalefactor(tau)/a_eq;
      R=0.75 * Omega_b/Omega_ph * a * a_eq;    
      dotalpha = y[1] - 3.0 * sqr(calHeq)/a/a* (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);
      fprintf(tcafile,"%d  %f  %e  %e  %e  %e  %e  %e\n",0,tau+DeltaTau,y[1],y[2],y[6],dotalpha,sigma,Pie);
    }
  
  tcaderivs(tau_handoff,y,dydt);    /*   call this to get proper value for sigma, Pie  */

  /*  store values of perturbations at first step (i=1, tau=160) of the grid  */
  hdot[1] = y[1];
  Delta0[1] = y[2];
  thetab[1]= y[6];
  a=scalefactor(tau_handoff)/a_eq;
  R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* the inverse of MB' and CLASS's  R   */
  alphadot[1] = hdot[1] - 3.0 * sqr(calHeq)/a/a* (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);

  Delta2[1]=sigma/2.0;
  Pi[1]=Pie;


    /* this initializes Delta2[2,3] and Pi[2,3] to match the 2nd-order TCA with the IE  */
  for(i=2;i<=3;i++) {
    odeint(y,9,taugrid[i-1],taugrid[i],1.0e-6,DeltaTau/10.0,0.0,&nok,&nbad,tcaderivs,rkqc);

    tcaderivs(taugrid[i],y,dydt);
    a=scalefactor(taugrid[i])/a_eq;
    R=0.75 * Omega_b/Omega_ph * a * a_eq;    
    dotalpha = y[1] - 3.0 * sqr(calHeq)/a/a* (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);

  
    fprintf(tcafile,"%d  %f  %e  %e  %e  %e  %e  %e\n",i,taugrid[i],y[1],y[2],y[6],dotalpha,sigma,Pie);

    hdot[i] = y[1];
    Delta0[i] = y[2];
    thetab[i]= y[6];
    alphadot[i] = dotalpha;
    Delta2[i]=sigma/2.0;
    Pi[i]=Pie;
  }


  /*****************************************************************************************/
  /*   Set up initial ansatz for Delta_2 and Pi  on the grid   */

  /*   read initial values from previous run    */

    xefile=fopen("iterative.in","r");
  
  while(fscanf(xefile,"%d  %lf %lf  %lf  %lf  %lf  %lf  %lf",
	       &i,&tau,&aa,&bb,&cc,&dd,&D2,&ppp) != EOF) {
    if((tau-taugrid[i]) > 0.001) {
	printf("tau= %f    but    taugrid[%d]=%f\n",tau,i,taugrid[i]);
	exit(1);
    }
    Delta2[i]=D2/2.0;
    Pi[i]=ppp;
  }
  fclose(xefile);
  

  /* Here is a pretty crappy analytic guesstimate ---
         does not seem to do any better than Delta2=Pi=0  
  if(taugrid[i]<=280.0) {
    Delta2[i] = 0.2*sin( r_s[i]*k)*sqr((taugrid[i]/280.0));
    rsmatch=r_s[i];
  }
  else {
    Delta2[i] = 0.2*sin( (rsmatch + sqrt(3.0)* (rsmatch/280.0)*(taugrid[i]-280.0)) * k)
      *exp(-sqr((taugrid[i]-280.0)/50.0));
    Pi[i]=Delta2[i];
  }
   */

  /*  Here we'll try setting Delta_2 and Pi to zero initially and see what happens  */
  /* 
     for(i=4;i<=imax;i++) {
     Delta2[i] = 0.0;
     Pi[i]=0.0;
     }
  */
/*      
	testfile=fopen("test.out","w");
	for(i=1;i<=imax;i++) {
	fprintf(testfile,"%f  %e  %e\n",taugrid[i],Delta2[i],r_s[i]);
	}
	fclose(testfile);
	exit(1);
*/


      
/*****************************************************************************************/
/*   now start the iteration, integrating DEs first and then IEs  */

for(iteration=1;iteration<=5;iteration++) {

  /*  first initialize the Delta2 and Pi interpolation tables */
  spline(taugrid,Delta2,Npts,1.0e40,1.0e40,Delta22_array);
  spline(taugrid,Pi,Npts,1.0e40,1.0e50,Pi2_array);
    /*
      testfile=fopen("test.out","w");
      for(tau=160.0;tau<=450.0;tau+=1.0) {
      fprintf(testfile,"%f  %f  %f\n",tau,Delta2_interp(tau),Pi_interp(tau));
      }
    */

      

    /*****************************************************************************************/
    /*  now integrate DEs, storing values of pert variables at taugrid as we go  */

    /*  first initialize perturbation variables to their values at end of TCA  */
    for(i=1;i<=9;i++) {
      ylate[i]=y[i];
    }

    for(i=4;i<=Npts;i++) {
      odeint(ylate,9,taugrid[i-1],taugrid[i],1.0e-6,0.1,0.0,&nok,&nbad,exactderivs,rkqc);

      /*  store values of perturbations needed for IE at each step of the grid  */
      hdot[i] = ylate[1];
      Delta0[i] = ylate[2];
      thetab[i]= ylate[6];
      a=scalefactor(taugrid[i])/a_eq;
      R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* the inverse of MB' and CLASS's  R   */
      alphadot[i] = hdot[i] - 3.0 * sqr(calHeq)/a/a* (Omega_ph/Omega_r)*((ylate[3] + R * ylate[6]) + fnu  * ylate[8]);
    }
    /* if(iteration==1) {
      testfile=fopen("test.out","w");
      for(i=1;i<=Npts;i++) {
	fprintf(testfile,"%d  %f  %e  %e  %e  %e\n",i,taugrid[i],hdot[i],Delta0[i],thetab[i],alphadot[i]);
      }
      fclose(testfile);
      exit(1);
    }
    */


    /*****************************************************************************************/
    /*  now solve integral equation */
    for(count=4;count<=Npts;count++) {    /*   evalute Delta2 and Pi at each time step */
	
      /*  initialize Delta2 and Pi   */
      Delta2[count]=0.0;
      Pi[count]=0.0;
      tau=taugrid[count];

      for(i=1;i<=count-1;i++) {     /* sum over all earlier times to evalute integrals */
	Deltakappa=(kappafrom1[i+1]-kappafrom1[i]);
	Deltakappa2=sqr(Deltakappa);

	/*  Kip1, Ki, weight contributions of integrands at i+1 and i  */
	if(Deltakappa>0.001) {
	  edk = exp(-Deltakappa);
	  Kip1=1.0 - edk - (1.0-(1.0+Deltakappa)*edk)/Deltakappa;
	  Ki = (1.0-(1.0+Deltakappa)*edk)/Deltakappa;
	}
	else {
	  Kip1 = Deltakappa/2.0;
	  Ki = Deltakappa/2.0;
	}

	Kip1 = Kip1 * exp(-kappafrom1[count]+kappafrom1[i+1]);
	Ki = Ki * exp(-kappafrom1[count]+kappafrom1[i+1]);

	/*  Iip1, Ii are the Delta2 integrands at i+1 and i
	    and Piip1, Pii are similar for Pi  */
	x = k * (tau-taugrid[i]);    /*********/
	Ii= ( -1.0/6.0*hdot[i]/kappadotprime[i] + Delta0[i]) * j2(x)
	  - (1.0/3.0 * alphadot[i]/kappadotprime[i]
	     +1.0/2.0*Pi[i]) * RLL2(x) + thetab[i]*j2prime(x)/k;

	Pii = 9.0*Pi[i]*j2overx2(x);
	/*   if i=count, we leave out the i+1 term here; it appears below given the semi-implicit nature of the integration  */

	x = k * (tau-taugrid[i+1]);       /*********/
	if(i<count-1) {
	  Piip1 = 9.0*Pi[i+1]*j2overx2(x);
	  Iip1 = ( -1.0/6.0*hdot[i+1]/kappadotprime[i+1] + Delta0[i+1]) * j2(x)
	    - (1.0/3.0 * alphadot[i+1]/kappadotprime[i+1]
	       +1.0/2.0*Pi[i+1]) * RLL2(x) + thetab[i+1]*j2prime(x)/k;
	}
	else {
	  Piip1=0.0;
	  Iip1=  ( -1.0/6.0*hdot[i+1]/kappadotprime[i+1] + Delta0[i+1]) * j2(x)
	    - (1.0/3.0 * alphadot[i+1]/kappadotprime[i+1]) * RLL2(x) + thetab[i+1]*j2prime(x)/k;
	  /*  Iip1 = (1.0/15.0)*alphadot[i+1]/kappadotprime[i+1];*/
	}
	Delta2[count] += Iip1*Kip1 + Ii*Ki;
	Pi[count] += Piip1*Kip1 + Pii*Ki;
      }      /*    end of sum over integrand  */

      /*   invert the 2x2 matrix from the semi-implicit integration scheme */
      sum = Pi[count]+Delta2[count];
      Pi[count] = sum /(1.0-(7.0/10.0)*Kip1);
      Delta2[count] += 0.1*Pi[count+1]*Kip1;

    }    /*  end of IEb   */

    if(iteration==1) {
      testfile=fopen("iterative1.out","w");
      for(i=4;i<=Npts;i++) {
	fprintf(testfile,"%d  %f  %e  %e  %e  %e  %e  %e\n",i,taugrid[i],hdot[i],Delta0[i],thetab[i],alphadot[i],2.0*Delta2[i],Pi[i]);
      }
      fclose(testfile);
    }

    if(iteration==2) {
      testfile=fopen("iterative2.out","w");
      for(i=4;i<=Npts;i++) {
	fprintf(testfile,"%d  %f  %e  %e  %e  %e  %e  %e\n",i,taugrid[i],hdot[i],Delta0[i],thetab[i],alphadot[i],2.0*Delta2[i],Pi[i]);
      }
      fclose(testfile);
    }

    if(iteration==3) {
      testfile=fopen("iterative3.out","w");
      for(i=4;i<=Npts;i++) {
	fprintf(testfile,"%d  %f  %e  %e  %e  %e  %e  %e\n",i,taugrid[i],hdot[i],Delta0[i],thetab[i],alphadot[i],2.0*Delta2[i],Pi[i]);
      }
      fclose(testfile);
    }

    if(iteration==4) {
      testfile=fopen("iterative4.out","w");
      for(i=4;i<=Npts;i++) {
	fprintf(testfile,"%d  %f  %e  %e  %e  %e  %e  %e\n",i,taugrid[i],hdot[i],Delta0[i],thetab[i],alphadot[i],2.0*Delta2[i],Pi[i]);
      }
      fclose(testfile);
    }

    if(iteration==5) {
      testfile=fopen("iterative5.out","w");
      for(i=4;i<=Npts;i++) {
	fprintf(testfile,"%d  %f  %e  %e  %e  %e  %e  %e\n",i,taugrid[i],hdot[i],Delta0[i],thetab[i],alphadot[i],2.0*Delta2[i],Pi[i]);
	fprintf(tcafile,"%d  %f  %e  %e  %e  %e  %e  %e\n",i,taugrid[i],hdot[i],Delta0[i],thetab[i],alphadot[i],2.0*Delta2[i],Pi[i]);
      }
      fclose(testfile);
    }
  }    /*  end of iteration  */
 fclose(tcafile);  

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


void tcaderivs(double t, double *y, double *dydt)
/*    differential equations for tight-coupling approximation  */
{
  double calH(),a,R,cH,dotkappa(),kappadot,thetabprime,slip,dotalpha,kappadd;
  double aprime_over_a,cHprime,thetadot0,alphadd,F,Fdot,sigmadot,slip1,sigma1;
  double ss2,cs2();
  
  /*   y[1] = hdot   y[2] = Delta0    y[3] =  thetag   y[4] = deltac  y[5] = deltaab  y[6]=thetab   */

  ss2=cs2(t);     /*    get the sound speed squared  */
  a=scalefactor(t)/a_eq;
  R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* note that this is the inverse of MB' and CLASS's  R   */
  cH = calH(t);
  cHprime = (calH(t*1.01)-calH(0.99* t))/(0.02*t);
  kappadot=dotkappa(t);
  
  /*  calculate these two first because they are needed for TCA quantities  */
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

    Pie=1.25*sigma - (1.0/3.0)*alphadd/kappadot - (2.0/3.0)*thetadot0/kappadot;      /*   this is the value of Pi in the TCA --- putting it here to feed to ICs for IE --- probably not very efficient   */

  
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
}


void exactderivs(double t, double *y, double *dydt)
/*    differential equations for tight-coupling approximation  */
{
  double calH(),a,R,cH,dotkappa(),kappadot,thetabprime;
  double Delta2_interp();
  double ss2,cs2();
  
  /*   y[1] = hdot   y[2] = Delta0    y[3] =  thetag   y[4] = deltac  y[5] = deltaab  y[6]=thetab   */

  ss2=cs2(t);     /*    get the sound speed squared  */
  a=scalefactor(t)/a_eq;
  R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* note that this is the inverse of MB' and CLASS's  R   */
  cH = calH(t);
  kappadot=dotkappa(t);

  dydt[1]= - cH*y[1] -3.0 * sqr(calHeq) * ( (4.0/a/a) * (Omega_ph/Omega_r) * ( y[2] + fnu*y[7]/4.0) + 0.5 * Omega_c/(Omega_b+Omega_c) /a * y[4] +0.5*y[5]*Omega_b/(Omega_b+Omega_c)/a);
  thetabprime = -cH*y[6]+ ss2*k*k*y[5] + kappadot/R*(y[3]-y[6]);
  dydt[2] = -y[3]/3.0 - y[1]/6.0;
  dydt[3] =  -R*(thetabprime+cH*y[6] - ss2*k*k*y[5])+ k*k*(y[2]-2.0*Delta2_interp(t));
  dydt[4] = -y[1]/2.0;
  dydt[5] = -y[6]-y[1]/2.0;
  dydt[6] = thetabprime;
  dydt[7] = -(4.0/3.0)*y[8] - (2.0/3.0)*y[1];
  dydt[8] = k*k*(y[7]/4.0 - y[9]);
  dydt[9] = -(3.0/t) * y[9] + (2.0/3.0)*y[8] + (1.0/3.0)* y[1];
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


  
double rsintegrand(double t)     /*   integrand for sound horizon --- used to initialize Delta_2  */
{
  double R;

  R=0.75 * Omega_b/Omega_ph  *scalefactor(t);    /* note that this is the inverse of MB' and CLASS's  R   */

  return sqrt( (1.0/3.0)/(1.0+R));
}


double Delta2_interp(double t)    /*   interpolated Delta2  */
{
  void splint();
  double y;

  splint(taugrid,Delta2,Delta22_array,Npts,t,&y);

  return y;
}


double Pi_interp(double t)    /*   interpolated Delta2  */
{
  void splint();
  double y;

  splint(taugrid,Pi,Pi2_array,Npts,t,&y);

  return y;
}


