#include<math.h>

void mydigamma(double * X, double * Y, double * coeffs);
void mytrigamma(double * X, double * Y, double * coeffs);
void myquadgamma(double * X, double * Y, double * coeffs);
void invtrigamma(double * X, double * Y, double * coeffs);
void sigmashrink(double * s2, double * d, double * s2t, double * d0, double * s20, double * gamcoeffs, int * n);
void getd0s20(double * s2, double * d, double * d0, double * e, double * s20, double * gamcoeffs, int * n);
void mymean(double * x, double * mean, int * n);


void sigmashrink(double * s2, double * d, double * s2t, double * d0, double * s20, double * gamcoeffs, int * n)
{
  int g;
  getd0s20(s2, d, d0, s2t, s20, gamcoeffs, n); // use s2t as buffer space for e
  if (*d0 > 0) for (g = 0; g < *n; g++) s2t[g] = (*d0*(*s20) + d[g]*s2[g])/(*d0 + d[g]);
  else for (g = 0; g < *n; g++) s2t[g] = *s20;
}

void getd0s20(double * s2, double * d, double * d0, double * e, double * s20, double * gamcoeffs, int * n)
{
  double X, Y;
  double meane;
  int g;
  for (g = 0; g < *n; g++)
  {
    X = d[g]/2;
    mydigamma(&X, &Y, gamcoeffs);
    e[g] = log(s2[g]) - Y + log(X);
  }  
  mymean(e, &meane, n);
  for (g = 0; g < *n; g++)
  {
    X = d[g]/2;
    mytrigamma(&X, &Y, gamcoeffs);
    e[g] = (e[g]-meane)*(e[g]-meane)*(*n)/(*n-1) - Y;
  }    
  mymean(e, &X, n); 
  invtrigamma(&X, &Y, gamcoeffs);
  *d0 = Y*2;
  if (*d0 > 0)
  {
    mydigamma(&Y, &X, gamcoeffs);
    *s20 = exp(meane + X - log(Y));
  }
  else *s20 = exp(meane);
}


void mymean(double * x, double * mean, int * n)
{
  int k;
  *mean = 0;
  for (k = 0; k < *n; k++) *mean = *mean + x[k];
  *mean = *mean / (*n);
  return;
}
























  void mydigamma(double * X, double * Y, double * coeffs)
  {
    double y;
    double x;
    double x2;
    double xpows;
    int i;
    x = *X;
    y = 0;
    while (x < 10)
    {
      y += 1/x;
      x += 1;  
    }
    x2 = x*x;
    xpows = 1/x2;
    y += coeffs[0]*xpows;
    for (i = 1; i < 20; i++) 
    {
      xpows = xpows/x2;
      y += coeffs[i]*xpows;
    }
    y = y + 1/(2*x);
    y = y - log(x);
    *Y = -y;
  }

  void mytrigamma(double * X, double * Y, double * coeffs)
  {
    double y, y0;
    double x;
    double x2;
    double xpows;
    int i;
    x = *X;
    y = y0 = 0;
    while (x < 10)
    {
      y0 += 1/(x*x);
      x += 1;  
    }
    x2 = x*x;
    xpows = 1/x2;
    y += coeffs[20]*xpows;
    for (i = 21; i < 40; i++) 
    {
      xpows = xpows/x2;
      y += coeffs[i]*xpows;
    }
    y = y / x;
    y = y + 1/(2*x2);
    y = y + 1/x;
    *Y = y + y0;
  }

  void myquadgamma(double * X, double * Y, double * coeffs)
  {
    double y, y0;
    double x;
    double x2;
    double xpows;
    int i;
    x = *X;
    y = y0 = 0;
    while (x < 10)
    {
      y0 += 1/(x*x*x);
      x += 1;  
    }
    x2 = x*x;
    xpows = 1/x2;
    y += coeffs[40]*xpows;
    for (i = 41; i < 60; i++) 
    {
      xpows = xpows/x2;
      y += coeffs[i]*xpows;
    }
    y = y / x2;
    y = y + 1/(2*x2*x);
    y = y + 0.5/x2;
    *Y = -2*(y+y0);
  }

  void invtrigamma(double * X, double * Y, double * coeffs)
  {
    double x;
    double y;
    double di;
    double tgy, qgy;
    x = *X;
    if (x > 10000000)
    {
      *Y = 1/sqrt(x);
    }
    else if (x < 0.000001)
    {
      *Y = 1/x;
    }
    else
    {
      y = 0.5 + 1/x;
      mytrigamma(&y, Y, coeffs);
      tgy = *Y;
      myquadgamma(&y, Y, coeffs);
      qgy = *Y;
      di = tgy*(1-tgy/x)/qgy;
      while (-di/y > 0.00000001)
      {
        y = y + di;
        mytrigamma(&y, Y, coeffs);
        tgy = *Y;
        myquadgamma(&y, Y, coeffs);
        qgy = *Y;
        di = tgy*(1-tgy/x)/qgy;
      }
      *Y = y;
    }
  }

