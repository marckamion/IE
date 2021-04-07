/*     time evolution of of perturbation variables with IE for one k mode  */

/*       now includes neutrino flucts using CLASS UFA   */
/*    TO DO:     
                 fix late-time integration  
		 add baryon sound speed
		 figure out allowable/optimal taustep
		 figure out allowable starting time
		 loop over k
		 compare with CLASS/CAMB output for pert variables after recombination  
                 figure out how to include HyRec?
                 figure out how to use exact evolution of scale factor with tau?  */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/Users/kamion/SCIENCE/recipes_c/nrutil.c"
#include "qromb.c"
#include "/Users/kamion/SCIENCE/recipes_c/trapzd.c"
#include "/Users/kamion/SCIENCE/recipes_c/polint.c"
#include "/Users/kamion/SCIENCE/recipes_c/splint.c"
#include "/Users/kamion/SCIENCE/recipes_c/spline.c"
#include "/Users/kamion/SCIENCE/recipes_c/bessj1.c"
#include "/Users/kamion/SCIENCE/recipes_c/odeint.c"
#include "/Users/kamion/SCIENCE/recipes_c/bsstep.c"
#include "/Users/kamion/SCIENCE/recipes_c/rkqc.c"
#include "/Users/kamion/SCIENCE/recipes_c/rk4.c"
#include "/Users/kamion/SCIENCE/recipes_c/mmid.c"
#include "/Users/kamion/SCIENCE/recipes_c/rzextr.c"

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
static double taustep,tau_handoff;
static double *z_array,*xe_array,*z2_array,*xe2_array,*Tm_array,*Tm2_array;

static double sigma;

int main(int argc,char *argv[])
{
  double tau0integrand(),xe(),scalefactor(),dotkappa(),calH();
  double j2overx2(),j2(),j2prime(),RLL2();
  double tau;   /*  conformal time in Moc   */
  double kappadot,kappa,visibility;
  double ktauinitial;    /* initial value of k*tau for TCA integration  */
  double tauinitial;
  double tau_end,Pie;
  double kt2,kt3,kt4,om,eps,Rb;
  double *y,*dydt;    /*    TCA perturbation variables  */
  double taumax;    /*  start and end of integration   */
  void odeint(),bsstep(),rkqc(),tcaderivs();
  int nok,nbad;
  int Npts;     /* number of integration time steps   */
  double *Delta0,*hdot,*alphadot,*thetag,*thetab,*deltab,*deltac,*Delta2,*Pi,*tauj,*kappafrom1,*kappadotprime;
  double R,a,x,Rold,cH,thetabprime,dotalpha,kappadd;
  double test1,test2,test3,dekappa,eip1,ei,Deltakappa1,Deltakappa;
  double Deltakappa2,Deltakappa4,Deltakappa3;
  double Kip1,Ki,Kim1,Iip1,Ii,Iim1,Pii,Piip1,Piim1,edk,dk,dk1,thetabpp,sum;
  double x0,x1,x2;
  double z,fraction,Tm;

  int i,j,count,n;

  FILE *tcafile,*kappadotfile,*testfile,*xefile;
  y=vector(1,9);      /*   number of TCA perturbation variables  */
  dydt=vector(1,9);      /*   number of TCA perturbation variables  */

  tcafile=fopen("tcaresults.out","w");
  kappadotfile=fopen("kapparesults.out","w");
  testfile=fopen("test.out","w");


  Npts=6000;
  Delta0=vector(1,Npts);
  Delta2=vector(1,Npts);
  hdot=vector(1,Npts);
  alphadot=vector(1,Npts);
  deltac=vector(1,Npts);
  deltab=vector(1,Npts);
  thetag=vector(1,Npts);
  thetab=vector(1,Npts);
  Pi=vector(1,Npts);
  tauj=vector(1,Npts);
  kappafrom1=vector(1,Npts);
  kappadotprime=vector(1,Npts);
  z_array=vector(1,8000);
  xe_array=vector(1,8000);
  z2_array=vector(1,8000);
  xe2_array=vector(1,8000);
  Tm_array=vector(1,8000);
  Tm2_array=vector(1,8000);



  /*   Planck best-fit values  */
  hhh=0.674;
  Omega_c=0.12038/hhh/hhh;
  Omega_b=0.0224/hhh/hhh;
  fnu=0.6904;
  Omega_ph=2.4735e-5/hhh/hhh;
  Omega_r=Omega_ph*(1.0+fnu);
  Hinverse=2998.0/hhh;
  T0=2.7255 * 8.6217e-5;
  rhocrit = 1.054e-5 * hhh * hhh;
  ne0 = rhocrit*Omega_b/0.938;    /*   this is actually the nucleon density  */

  a_eq = Omega_r/(Omega_c + Omega_b);

  /* calculate conformal time today  */
  tau0=qromb(tau0integrand,0.0000000000001,1.0) * Hinverse;
  tau_eq=qromb(tau0integrand,0.0000000000001,Omega_r/(Omega_c+Omega_b)) * Hinverse;
  tau_ls=qromb(tau0integrand,0.0000000000001,1.0/1100.0) * Hinverse;
  printf("tau0=%f     tau_eq=%f  tau_ls=%f  z_eq=%f \n",tau0,tau_eq,tau_ls,1.0/(Omega_r/(Omega_c+Omega_b)));
  calHeq=calH(tau_eq);

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
  kappa=0.0;
  for(tau = tau_ls*5.0; tau>=tau_ls/20.0; tau -= (5.0-1.0/20.0)*tau_ls/1000.) {
    kappadot=dotkappa(tau);
    kappa += kappadot * 4.5*tau_ls/1000.0;
    visibility = kappadot*exp(-kappa);
    fprintf(kappadotfile,"%f   %f   %f   %f   %f\n",tau,1.0/scalefactor(tau),xe(tau), 1.0/kappadot,visibility); 
  }
  fclose(kappadotfile);


  for(k=1.0/tau0; k<= 3001.00001/tau0; k+=30.0/tau0) {
      
    /*    k = 20.0/200.0;  */   /*   wavenumber  */
  

    /*   integrate tight-coupling eqns from tauinitial to tau_handoff  */

    tauinitial=10.0;
    tau_handoff=160.0;
    tau_end = 450.0;
    taustep = 0.25;      /*  step forward in conformal-time steps of 1.0 Mpc  */


    /*   y[1] = hdot   y[2] = Delta0    y[3] =  thetag   y[4] = deltac  y[5] = deltaab  y[6]=thetab     y[7] = delta_nu   y[8]=theta_nu   y[9]=sigma_nu   using CLASS eqns   */

    /*   initial conditions   from Cyr-Racine & Sigurdson   */

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
    /*  fprintf(tcafile,"%f  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e\n",tauinitial,1.0/(kappadot*(tauinitial+taustep)),kappadot*taustep,y[1],y[2],y[3],y[4],y[5],y[6],0.0,0.0);   */        /*   print initial conditions   */

    count=0;  kappa=0.0;
    for(tau=tauinitial; tau<=tau_end; tau+=taustep) {
      if(tau < tau_handoff) {
	odeint(y,9,tau,tau+taustep,1.0e-6,taustep/10.0,0.0,&nok,&nbad,tcaderivs,rkqc);

	tcaderivs(tau+taustep,y,dydt);    /*   call this to get proper value for sigma  */

	/*   start collecting initial data for IE  */
	if(tau > tau_handoff - 20.0) {
	  count++;
	  kappadot=dotkappa(tau);
	  kappa+=kappadot*taustep;
	  Delta0[count]=y[2];
	  Delta2[count]=sigma/2.0;
	  hdot[count]=y[1];
	  deltac[count]=y[4];
	  deltab[count]=y[5];
	  thetag[count]=y[3];
	  thetab[count]=y[6];
	  Pi[count]= 1.25*sigma;        /*  need to check signs on Pi   */
	  tauj[count]=tau+taustep;
	  kappafrom1[count]=kappa;
	  kappadotprime[count]=kappadot;
	  a=scalefactor(tau)/a_eq;
	  R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* note that this is the inverse of MB's and CLASS's  R   */
	  /*	  alphadot[count] = hdot[count] - 3.0 * sqr(calHeq)/a/a*(thetag[count] + R * thetab[count]);     without neutrinos  */
	  alphadot[count] = hdot[count] - 3.0 * sqr(calHeq)/a/a* (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);   /*  with neutrinos  */
	}
	Pie=1.25*sigma;
      }

      /*  do "exact" integral equations  */
      else {
	/*   calculate the derivatives and advance solutions of ODEs  */
	sigma = 2.0 * Delta2[count];
	tcaderivs(tau,y,dydt);
	for(i=1;i<=9;i++) {
	  y[i] += dydt[i]*taustep;
	}
	a=scalefactor(tau)/a_eq;
	R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* note that this is the inverse of MB's and CLASS's  R   */
	kappadot=dotkappa(tau);
	kappa+=kappadot*taustep;
  
	hdot[count+1] = y[1];
	Delta0[count+1] = y[2];
	thetag[count+1] = y[3];
	thetab[count+1]= y[6];
	deltac[count+1] = y[4];
	deltab[count+1] = y[5];
	tauj[count+1]=tau+taustep;
	kappafrom1[count+1]=kappa;
	kappadotprime[count+1]=kappadot;
	/*   alphadot[count+1] = hdot[count+1] - 3.0 * sqr(calHeq)/a/a*(thetag[count+1] + R * thetab[count+1]);   without neutrinos   */
	  alphadot[count+1] = hdot[count+1] - 3.0 * sqr(calHeq)/a/a* (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);   /*  with neutrinos  */
	/*    now do the integral equations   */
	/*  initialize Delta2 and Pi   */
	Delta2[count+1]=0.0;
	Pi[count+1]=0.0;

	/*    test1 = -exp(-kappafrom1[count+1]+kappafrom1[1]);  */
	test1 = 0.0;
	test2 = 0.0;
	test3 = 0.0;
	for(i=2;i<=count;i++) {
	  x=k * (tau-tauj[i]);
	  Deltakappa=(kappafrom1[i+1]-kappafrom1[i]);
	  Deltakappa2=Deltakappa*Deltakappa;
	  Deltakappa3=Deltakappa2*Deltakappa;
	  Deltakappa4=Deltakappa2*Deltakappa2;

	  if(Deltakappa>0.01) {
	      
	    edk = exp(-Deltakappa);
	    Kip1=1.0 - edk - (1.0-(1.0+Deltakappa)*edk)/Deltakappa + (1.0 - (1.0+Deltakappa + (0.5)*Deltakappa2)*edk)/Deltakappa2;
	    Ki = (1.0-(1.0+Deltakappa)*edk)/Deltakappa - 2.0* (1.0 -(1.0+Deltakappa + 0.5*Deltakappa2)*edk)/Deltakappa2;
	    Kim1 =  (1.0-(1.0+Deltakappa+0.5*Deltakappa2)*edk)/Deltakappa2;
	  }
	  else {
	    Kip1 = 2.0*Deltakappa/3.0 -7.0*Deltakappa2/24.0 + 11.0*Deltakappa3/120.0;
	    Ki = Deltakappa/6.0 - Deltakappa2/12.0 + Deltakappa3/40.0;
	    Kim1 = Deltakappa/6.0 - Deltakappa2/8.0 + Deltakappa3/20.0;
	    }

	  /*   x2=kappafrom1[i+1];  x1=kappafrom1[i];  x0=kappafrom1[i-1];
	       dk=x2-x1;  dk1=x1-x0;  */
	  
	  /*  Kip1=((x1-x0-2.0)*edk+2.0+x1+x0*(1.0+x1-x2)-2.0*x2-x1*x2+x2*x2)/(x2-x0)/(x2-x1);
	      Ki = (x2-x0-2.0+edk*(2.0+x1*x1+x2-x1*(2.0+x2)+x0*(1.0-x1+x2)))/(x0-x1)/(x1-x2);
	      Kim1 = (edk*(-2.0+x1-x2)+2.0+x1-x2)/(x0-x1)/(x0-x2);  */

	  /*  Kip1=((dk1-2.0)*edk+2.0+x1+x0*(1.0-dk)-2.0*x2-x1*x2+x2*x2)/(x2-x0)/dk;
	      Ki = (x2-x0-2.0+edk*(2.0+x1*x1+x2-x1*(2.0+x2)+x0*(1.0-dk)))/dk1/dk;
	      Kim1 = (edk*(-2.0-dk)+2.0-dk)/dk1/(dk+dk1);   */

	  Kip1 = Kip1 * exp(-kappafrom1[count+1]+kappafrom1[i+1]);
	  Ki = Ki * exp(-kappafrom1[count+1]+kappafrom1[i+1]);
	  Kim1 = Kim1 * exp(-kappafrom1[count+1]+kappafrom1[i+1]);

	  Iim1= ( -1.0/6.0*hdot[i-1]/kappadotprime[i-1] + Delta0[i-1]) * j2(x+k*taustep)
	    - (1.0/3.0 * alphadot[i-1]/kappadotprime[i-1]
	       +1.0/2.0*Pi[i-1]) * RLL2(x+k*taustep) + thetab[i-1]*j2prime(x+k*taustep)/k;
	  Ii= ( -1.0/6.0*hdot[i]/kappadotprime[i] + Delta0[i]) * j2(x)
	    - (1.0/3.0 * alphadot[i]/kappadotprime[i]
	       +1.0/2.0*Pi[i]) * RLL2(x) + thetab[i]*j2prime(x)/k;

	  Piim1 = 9.0*Pi[i-1]*j2overx2(x+k*taustep);
	  Pii = 9.0*Pi[i]*j2overx2(x);
	  if(i<count) {
	    Piip1 = 9.0*Pi[i+1]*j2overx2(x-k*taustep);
	    Iip1 = ( -1.0/6.0*hdot[i+1]/kappadotprime[i+1] + Delta0[i+1]) * j2(x-k*taustep)
	      - (1.0/3.0 * alphadot[i+1]/kappadotprime[i+1]
		 +1.0/2.0*Pi[i+1]) * RLL2(x-k*taustep) + thetab[i+1]*j2prime(x-k*taustep)/k;
	  }
	  else {
	    Piip1=0.0;
	    Iip1=  ( -1.0/6.0*hdot[i+1]/kappadotprime[i+1] + Delta0[i+1]) * j2(x-k*taustep)
	      - (1.0/3.0 * alphadot[i+1]/kappadotprime[i+1]) * RLL2(x-k*taustep) + thetab[i+1]*j2prime(x-k*taustep)/k;
	    /*  Piip1=Pii;
		Iip1=( -1.0/6.0*hdot[i+1]/kappadotprime[i+1] + Delta0[i+1]) * j2(x-k*taustep)
		- (1.0/3.0 * alphadot[i+1]/kappadotprime[i+1]+0.5*Pi[i]) * RLL2(x-k*taustep) + thetab[i+1]*j2prime(x-k*taustep)/k;  */

	  }

	  Delta2[count+1] += Iip1*Kip1 + Ii*Ki + Iim1*Kim1;
	  Pi[count+1] += Piip1*Kip1 + Pii*Ki + Piim1*Kim1;
	  /*	  test1 += Kip1+Ki+Kim1;
	  test2 += -(Kip1*(kappafrom1[i+1]-kappafrom1[count+1]) + Ki*(kappafrom1[i]-kappafrom1[count+1]) + Kim1*(kappafrom1[i-1]-kappafrom1[count+1]));
	  test3 += 0.5*(Kip1*sqr(kappafrom1[i+1]-kappafrom1[count+1]) + Ki*sqr(kappafrom1[i]-kappafrom1[count+1]) + Kim1*sqr(kappafrom1[i-1]-kappafrom1[count+1]));
	  */

	  /*  printf("%d  %d  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",count,i,kappafrom1[i+1],kappafrom1[i],Deltakappa,Kip1,Ki,Kim1,Iip1,Ii,Iim1,Piip1,Pii,Piim1);   */
	  /*    printf("%d  %f  %f\n,",i,kappafrom1[i+1],kappafrom1[i]);  */

	}     /*    end of sum over integrand  */


	/*    Pi[count+1] += Delta2[count+1];   */
	sum = Pi[count+1]+Delta2[count+1];
	Pi[count+1] = sum /(1.0-(7.0/10.0)*Kip1);
	Delta2[count+1] += 0.1*Pi[count+1]*Kip1;

	sigma = 2.0*Delta2[count+1];
	Pie=Pi[count+1];
	count++;
      }     /*   end of late-time evaluation of integro-differential eqns  */

      fprintf(tcafile,"%1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.3e  %1.6e  %1.6e  %1.3e  %1.3e\n",k,tau+taustep,y[1],y[2],y[3],y[4],y[5],y[6],sigma,Pie,y[3]-y[6],
	      (1.0/k/k/k)*(Omega_c *y[4]+Omega_b*y[5])*(Omega_c *y[4]+Omega_b*y[5]) );

      fprintf(testfile,"%1.3e  %d  %f  %f  %f  %f\n",k,count,tauj[count],alphadot[count],kappafrom1[count],kappadotprime[count]);

      printf("k=%f   tau=%f\n",k,tau);
      
      /*  printf("%f  %f  %f  %f\n",tau+taustep,test1,test2,test3); */
    }     /*   end of loop over tau from tauinitial to taufinal  */
  }
}




      
  
/*

  double xe(double tau)  
  {
  double a,A,scalefactor();
 
  return 0.5*(tanh((260.0-tau)/50.0)+1.0);    
  }

*/


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



  
double conformaltime(double a)
{
  double sqrt2m1=0.4142135623731;

  return tau_eq * (sqrt( (a/a_eq)+1.0) - 1.0)/sqrt2m1;
}
  

double dotkappa(double tau)    /*   \dot\kappa in units of 1/Mpc  */
{
  double xe(),a,scalefactor();
  static double YHe;         /*   helium fraction   */

  a=scalefactor(tau);

  return 2.052 * ne0 * (1.0-YHe) * xe(tau) / a/a;
}
  


  
double tau0integrand(double a)
{
  double E();

  return 1.0/a/a/E(1.0/a-1.0);
}


double j1(double x)
{
  return sin(x)/x/x - cos(x)/x;
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
 
  /*   y[1] = hdot   y[2] = Delta0    y[3] =  thetag   y[4] = deltac  y[5] = deltaab  y[6]=thetab   */

  a=scalefactor(t)/a_eq;
  R=0.75 * Omega_b/Omega_ph * a * a_eq;    /* note that this is the inverse of MB' and CLASS's  R   */
  cH = calH(t);
  cHprime = (calH(t*1.01)-calH(0.99* t))/(0.02*t);
  kappadot=dotkappa(t);
  
  /*  calculate these two first because they are needed for TCA quantities  */
  dydt[1]= - cH*y[1] -3.0 * sqr(calHeq) * ( (4.0/a/a) * (Omega_ph/Omega_r) * ( y[2] + fnu*y[7]/4.0) + 0.5 * Omega_c/(Omega_b+Omega_c) /a * y[4] +0.5*y[5]*Omega_b/(Omega_b+Omega_c)/a);
  dydt[8] = k*k*(y[7]/4.0 - y[9]);

  /*   use TCA for sigma and slip at early times  */
  if(t<tau_handoff+taustep) {
    dotalpha = y[1] - 3.0 * sqr(calHeq)/a/a * (Omega_ph/Omega_r)*((y[3] + R * y[6]) + fnu  * y[8]);
    
    aprime_over_a = cH*cH+cHprime;
    kappadd = (dotkappa(1.01*t)-dotkappa(0.99*t))/(0.02*t);
    thetadot0 = R/(1.0+R) *(-cH*y[3] + k*k/R*y[2]);

    /*   alphadd = dydt[1] - 3.0*sqr(calHeq) * ( y[3]*cH*(R/(1.0+R)-2.0)+thetadot0);  this is what I had before, but I don't think it looks right   */
    alphadd = dydt[1] - 2.0*cH * (dotalpha-y[1]) - 3.0*sqr(calHeq)/a/a * Omega_ph/Omega_r * (thetadot0*(1.0+R) + cH*R*y[3] + fnu * dydt[8]);
    
  
    /*  lowest-order TCA approx for sigma (photon quadrupole) from CLASS paper  */
    sigma = 8.0/45.0/kappadot*(2.0 * y[3]+ dotalpha);
    sigma1= sigma;

    /*  now get sigma to second order in 1/kappadot  */
    sigma = (1.0+11.0/6.0*kappadd/kappadot/kappadot)*sigma - (11.0/6.0)/kappadot *  (8.0/45.0)/kappadot* (2.0*thetadot0+alphadd);
  
    /* lowest-order TCA approx for slip  */
    F = R/(1.0+R)/kappadot;
    slip = -(kappadd/kappadot + 2.0*cH*R/(1.0+R))*(y[3]-y[6]) - F * (-aprime_over_a * y[6] + k*k*(-2.0*cH*y[2] - ( -y[3]/3.0-y[1]/6.0)));
    slip1 = slip;

    /*  now get slip to partial 2nd order (compromise_class approx)  */
    sigmadot = -kappadd/kappadot*sigma+(8.0/45.0)/kappadot*(2.0*thetadot0+alphadd);
    Fdot= 1.0/(1.0+R)/(1.0+R)*cH*R -R/(1.0+R)*kappadd/kappadot/kappadot;
    slip = (1.0- 2.0*cH*F)*slip - F*k*k*(2.0*cH*sigma+sigmadot-(1.0/3.0)*(F*thetadot0 +  2.0* Fdot*y[6]));
  }
  
  dydt[2] = -y[3]/3.0 - y[1]/6.0;

  if(kappadot/R>10.0*cH && t<tau_handoff+taustep) {
    /* TCA calculation of thetabprime   */
    F = R/(1.0+R)/kappadot;
    slip = -(kappadd/kappadot + 2.0*cH*R/(1.0+R))*(y[3]-y[6]) - F * (-aprime_over_a * y[6] + k*k*(-2.0*cH*y[2] - ( -y[3]/3.0-y[1]/6.0)));
    slip1 = slip;

    /*  now get slip to partial 2nd order (compromise_class approx)  */
    sigmadot = -kappadd/kappadot*sigma+(8.0/45.0)/kappadot*(2.0*thetadot0+alphadd);
    Fdot= 1.0/(1.0+R)/(1.0+R)*cH*R -R/(1.0+R)*kappadd/kappadot/kappadot;
    slip = (1.0- 2.0*cH*F)*slip - F*k*k*(2.0*cH*sigma+sigmadot-(1.0/3.0)*(F*thetadot0 +  2.0* Fdot*y[6]));
    thetabprime= - R /(1.0+R) * (cH*y[6]-k*k/R*(y[2]-sigma)+slip/R);
  }
  else {
    /*  exact thetabprime  */
    thetabprime = -cH*y[6]+kappadot/R*(y[3]-y[6]);
  }
  dydt[3] =  -R*(thetabprime+cH*y[6])+k*k*(y[2]-sigma);
  dydt[4] = -y[1]/2.0;
  dydt[5] = -y[6]-y[1]/2.0;
  dydt[6] = thetabprime;
  dydt[7] = -(4.0/3.0)*y[8] - (2.0/3.0)*y[1];
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


  
