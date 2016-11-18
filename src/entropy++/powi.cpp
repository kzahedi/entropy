#include "powi.h"

int powi(int a, int b)
{
  if(b == 0) return 1;
  int r = 1;
  for(int i = 0; i < b; i++)
  {
    r *= a;
  }
  return r;
}

