#include "powi.h"

int powi(int a, int b)
{
  int r = 1;
  for(int i = 0; i < b; i++)
  {
    r *= a;
  }
  return r;
}

