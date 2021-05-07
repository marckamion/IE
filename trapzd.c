#define FUNC(x) ((*func)(x))

double trapzd(func,a,b,n)
double a,b;
double (*func)();                           /* ANSI: double (*func)(double); */
int n;
{
   double x,tnm,sum,del;
   static double s;
   static int it;
   int j;

   if (n == 1) {
      it=1;
      return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
   } else {
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
      it *= 2;
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}
