#include <entropy++/continuous/Functions.h>

#include <entropy++/defs.h>

#include <math.h>

double entropy::continuous::dist(vector<double> &a, vector<double> &b, vector<int> &indices)
{
  double d = 0.0;
  double v = 0.0;
  for(vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
  {
    v = a[*i] - b[*i];
    d += v * v;
  }
  return sqrt(d);
}

double entropy::continuous::maxNorm3(vector<double>& a, vector<double>& b, vector<int>& x, vector<int>& y, vector<int>& z)
{
  double xDist = dist(a, b, x);
  double yDist = dist(a, b, y);
  double zDist = dist(a, b, z);
  return MAX3(xDist, yDist, zDist);
}

double entropy::continuous::maxNorm2(vector<double>& a, vector<double>& b, vector<int>& x, vector<int>& y)
{
  double xDist = dist(a, b, x);
  double yDist = dist(a, b, y);
  return MAX(xDist, yDist);
}

double entropy::continuous::getEpsilon(int rowIndex,
                  entropy::DContainer *data,
                  vector<int> &xIndices,
                  vector<int> &yIndices,
                  vector<int> &zIndices,
                  int k)
{
  vector<double> distances(data->rows(), 0.0);
  vector<double> reference = data->row(rowIndex);
  vector<double> vrow;
  for(int row = 0; row < data->rows(); row++)
  {
    vrow = data->row(row);
    distances[row] = maxNorm3(reference, vrow, xIndices, yIndices, zIndices);
  }
  sort(distances.begin(), distances.end());
  return distances[k];
}

double entropy::continuous::harmonic(int n)
{
  double r = -0.5772156649;
  for(int i = 1; i <= n; i++)
  {
    r -= 1.0 / (float)i;
  }
  return r;
}

int entropy::continuous::count2(double epsilon,
           vector<double>& xyz,
           entropy::DContainer* data,
           vector<int>& xIndices,
           vector<int>& yIndices)
{
  int c = -1; // because we will also count xyz[t] vs. xyz[t]

  for(int t = 0; t < data->rows(); t++)
  {
    vector<double> vrow = data->row(t);
    if(maxNorm2(xyz, vrow, xIndices, yIndices) < epsilon) c++;
  }

  return c;
}

int entropy::continuous::count1(double epsilon,
           vector<double>& xyz,
           entropy::DContainer* data,
           vector<int>& indices)
{
  int c = -1; // because we will also count xyz[t] vs. xyz[t]

  for(int t = 0; t < data->rows(); t++)
  {
    vector<double> vrow = data->row(t);
    if(dist(xyz, vrow, indices) < epsilon) c++;
  }

  return c;
}


