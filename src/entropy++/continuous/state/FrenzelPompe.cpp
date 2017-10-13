#include "FrenzelPompe.h"

#include <entropy++/Container.h>
#include <entropy++/continuous/Functions.h>

#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;

double entropy::continuous::state::FrenzelPompe(entropy::DContainer* xyz,
                                                vector<int> xIndices,
                                                vector<int> yIndices,
                                                vector<int> zIndices,
                                                int k)
{
  double r = 0.0;

  double hk = harmonic(k - 1);

  for(int t = 0; t < xyz->rows(); t++)
  {
    vector<double> vrow = xyz->row(t);
    double epsilon = getEpsilon(t, xyz, xIndices, yIndices, zIndices, k);

    int cNxy = count2(epsilon, vrow, xyz, xIndices, yIndices);
    double hNxy = harmonic(cNxy);

    int cNyz = count2(epsilon, vrow, xyz, yIndices, zIndices);
    double hNyz = harmonic(cNyz);

    int cNz = count1(epsilon, vrow, xyz, zIndices);
    double hNz = harmonic(cNz);

    r += hNxy + hNyz - hNz;
  }

  r /= float(xyz->rows());

  r -= hk;

  return r;
}
