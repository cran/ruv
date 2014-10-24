#include<math.h>

void quicksort(double * data, int a, int b);
void bsweep(double * A, int * ind, int * p, int * N);
void getsizeratios(double * A, double * B, double * Aj, double * Bj, double * bycx, double * ac, double * ajtB, double * sizeratios, double * scaledbetac, double * bwx, int * K1, int * nc);

void quicksort(double * data, int a, int b)
{
  int i, j;
  double x, y;
  if (a < b) 
  {
    x = data[b];
    i = a - 1;
    for (j = a; j < b; j++)
    {
      if (data[j] <= x)
      {
        i++;
        y = data[i];
        data[i] = data[j];
        data[j] = y;
      }
    }
    y = data[i+1];
    data[i+1] = data[b];
    data[b] = y;
    quicksort(data, a, i);
    quicksort(data, i+2, b);
  }
} 

void bsweep(double * A, int * ind, int * p, int * N)
{
  int i,j,k,count;
  double d, b;
  for(count = 0; count < *N; count++)
  {
    k = ind[count];  
    d = A[k*(*p)+k]; 
    for (i = 0; i<*p; i++) A[i*(*p)+k] /= d;  
    for (i = 0; i<k; i++) 
    {
      b = A[k*(*p)+i]; 
      for (j = 0; j < *p; j++) A[j*(*p)+i] -= b*A[j*(*p)+k]; 
      A[k*(*p)+i] = -b/d; 
    }
    for (i = k+1; i<*p; i++)
    {
      b = A[k*(*p)+i];
      for (j = 0; j < *p; j++) A[j*(*p)+i] -= b*A[j*(*p)+k];
      A[k*(*p)+i] = -b/d;
    }
    A[k*(*p)+k] = 1/d;
  }
}


void getsizeratios(double * A, double * B, double * Aj, double * Bj, double * bycx, double * ac, double * ajtB, double * sizeratios, double * scaledbetac, double * bwx, int * K1, int * nc)
{
  int i, j, k, l;
  double betachat, divisor, denominator;
  for (i = 1; i <= *K1; i++)
  {
    k=i-1;
    l=1;
    bsweep(B, &k, K1, &l);
    for (j = 0; j < *nc; j++)
    {
      for (k = 0; k < i; k++) Aj[k] = A[k] - bycx[j]*ac[*K1*j + k];
      for (k = 0; k < i; k++)
      {
        ajtB[k] = 0;
        for (l = 0; l < i; l++) ajtB[k] = ajtB[k] + ac[*K1*j + l]*B[*K1*k + l];
      }
      denominator = 1;
      for (k = 0; k < i; k++) denominator = denominator - ajtB[k]*ac[*K1*j + k];
      for (k = 0; k < i; k++)
      {
        for (l = 0; l < i; l++) Bj[i*k + l] = B[*K1*k + l] + ajtB[k]*ajtB[l]/denominator;
      }
      for (k = 0; k < i; k++)
      {
        bwx[k] = 0;
        for (l = 0; l < i; l++) bwx[k] = bwx[k] + Aj[l]*Bj[i*k + l];
      }
      betachat = bycx[j];
      for (k = 0; k < i; k++) betachat = betachat - bwx[k]*ac[*K1*j + k];
      divisor = 0;
      for (k = 0; k < i; k++) divisor = divisor + bwx[k]*bwx[k];
      divisor = sqrt(1 + divisor);
      scaledbetac[j] = betachat / divisor;     
      scaledbetac[j] = ac[*K1*j + i-1]/scaledbetac[j]; ////// single sizeratio 
      if (scaledbetac[j] < 0) scaledbetac[j] = -scaledbetac[j];
    }
    quicksort(scaledbetac, 0, *nc-1);
    if (*nc % 2 == 1) sizeratios[i-1] = scaledbetac[*nc/2];
    else sizeratios[i-1] = (scaledbetac[*nc/2-1] + scaledbetac[*nc/2])/2; 
  }
}


