#include <math.h>

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP JMAX+1
#define K 5

double qromb(double (*func)(double),double a,double b)
{
   double ss,dss,trapzd();
   double s[JMAXP+1],h[JMAXP+1];
   int j;
   void polint(),nrerror();

   h[1]=1.0;
   for (j=1;j<=JMAX;j++) {
      s[j]=trapzd(func,a,b,j);
      if (j >= K) {
         polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
         if (fabs(dss) < EPS*fabs(ss)) return ss;
      }
      s[j+1]=s[j];
      h[j+1]=0.25*h[j];
   }
   nrerror("Too many steps in routine QROMB");  return 0;
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K
