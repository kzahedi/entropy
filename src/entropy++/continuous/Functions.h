#ifndef __CONTINUOUS_DEFS_H__
#define __CONTINUOUS_DEFS_H__

#include <entropy++/Container.h>

#include <vector>

using namespace std;

namespace entropy
{
  namespace continuous
  {
    double dist(vector<double> &a, vector<double> &b, vector<int> &indices);

    double maxNorm3(vector<double>& a, vector<double>& b, vector<int>& x, vector<int>& y, vector<int>& z);

    double maxNorm2(vector<double>& a, vector<double>& b, vector<int>& x, vector<int>& y);

    double getEpsilon(int rowIndex,
                      entropy::DContainer *data,
                      vector<int> &xIndices,
                      vector<int> &yIndices,
                      vector<int> &zIndices,
                      int k);

    double harmonic(int n);

    int count2(double epsilon,
               vector<double>& xyz,
               entropy::DContainer* data,
               vector<int>& xIndices,
               vector<int>& yIndices);

    int count1(double epsilon,
               vector<double>& xyz,
               entropy::DContainer* data,
               vector<int>& indices);
  }
}

#endif // __CONTINUOUS_DEFS_H__
