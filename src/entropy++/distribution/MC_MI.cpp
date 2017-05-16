#include "MC_MI.h"

#include <entropy++/distribution/MI.h>

double entropy::distribution::MC_MI(double**** pw2w1s1a1, int dimW2, int dimW1, int dimS1, int dimA1)
{
  double** pw2w1 = new double*[dimW2];
  for(int i = 0; i < dimW2; i++)
  {
    pw2w1[i] = new double[dimW1];
    for(int j = 0; j < dimW1; j++)
    {
      pw2w1[i][j] = 0.0;
    }
  }

  double** pa1s1 = new double*[dimA1];
  for(int i = 0; i < dimA1; i++)
  {
    pa1s1[i] = new double[dimS1];
    for(int j = 0; j < dimS1; j++)
    {
      pa1s1[i][j] = 0.0;
    }
  }

  for(int w2 = 0; w2 < dimW2; w2++)
  {
    for(int w1 = 0; w1 < dimW1; w1++)
    {
      for(int s1 = 0; s1 < dimS1; s1++)
      {
        for(int a1 = 0; a1 < dimA1; a1++)
        {
          pw2w1[w2][w1] += pw2w1s1a1[w2][w1][s1][a1];
          pa1s1[a1][s1] += pw2w1s1a1[w2][w1][s1][a1];
        }
      }
    }
  }

  return entropy::distribution::MI(pw2w1, dimW2, dimW1) -
         entropy::distribution::MI(pa1s1, dimA1, dimS1);
}

