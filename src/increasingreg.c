#include <stdlib.h>

void increasingreg(double *y, int *n);
void collapse(int *start, int *stop, double *sum, double *avg, int *n0);

void increasingreg(double *y, int *n)
{
  int i, n0;
  int *start, *stop;
  double *sum;
  double *avg;
  start = (int *) malloc(*n*sizeof(int));
  stop = (int *) malloc(*n*sizeof(int));
  sum = (double *) malloc(*n*sizeof(double));
  avg = (double *) malloc(*n*sizeof(double));
  n0 = 0;
  for (i = 0; i < *n; i++)
  {
    start[n0] = i;
    stop[n0] = i;
    sum[n0] = y[i];
    avg[n0] = y[i];
    n0++;
    collapse(start, stop, sum, avg, &n0);
  }
  n0 = 0;
  for (i = 0; i < *n; i++)
  {
    y[i] = avg[n0];
    if (i == stop[n0]) n0++;
  }
  free(start);
  free(stop);
  free(sum);
  free(avg);
  return;
}

void collapse(int *start, int *stop, double *sum, double *avg, int *n0)
{
  if (*n0 == 1) return;
  else if (avg[*n0 - 1] < avg[*n0 - 2])
  {
    sum[*n0-2] = sum[*n0-2] + sum[*n0-1];
    avg[*n0-2] = sum[*n0-2] / (stop[*n0-1] - start[*n0-2] + 1);
    stop[*n0-2] = stop[*n0-1];
    *n0 = *n0 - 1;
    collapse(start, stop, sum, avg, n0);
    return;
  }
  else return;
}
