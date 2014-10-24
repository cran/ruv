#include<stdio.h>

void wipeout(double *x, int *indicator, int *n, int *bin)
{
  int start, i, maxindex;
  double max;
  for (i = 0; i < *n; i++) indicator[i] = 1;
  start = 0;
  while (*n - start > *bin)
  {
    max = 0;
    for (i = 0; i < *bin; i++)
    {
      if (x[start + i] > max)
      {
        max = x[start + i];
        maxindex = start + i;
      }
    }
    indicator[maxindex] = 0;
    start = start + *bin;
  }
  max = 0;
  i = start;
  while (i < *n)
  {
    if (x[i] > max)
    {
      max = x[i];
      maxindex = i;
    }
    i++;
  }
  indicator[maxindex] = 0;
  return;
}
