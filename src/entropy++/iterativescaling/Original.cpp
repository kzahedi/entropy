#include "Original.h"
#include <entropy++/powi.h>

#include <iostream>

using namespace entropy::iterativescaling;
using namespace std;

Original::Original()
{
  _sizeAlphabet = 0;
  _alphabet     = NULL;
}

Original::Original(int n, vector<vector<int> > features, vector<double> p)
{
  _sizeAlphabet = pow(2,n);
  _alphabet     = new Matrix(_sizeAlphabet,n);
  __generateAlphabet(n);

  assert(p.size() == _sizeAlphabet);
  _targetp  = p;

  _features = features;

  int max = 0;
  for(int i = 0; i < (int)_features.size(); i++)
  {
    if(_features[i].size() > max) max = _features.size();
  }

  cout << "Expected number of iterations: " << max * _sizeAlphabet << endl;

  _p1 = vector<double>(_sizeAlphabet);
  _p2 = vector<double>(_sizeAlphabet);

  double val = 1.0 / _sizeAlphabet;
  for(int i = 0; i < _sizeAlphabet; i++)
  {
    _p1[i] = val;
  }
}

Original::~Original()
{
  delete _alphabet;
}

void Original::__generateAlphabet(int n)
{
  for(int z = 0; z < _sizeAlphabet; z++)
  {
    for(int s = n; s > 0; s--)
    {
      (*_alphabet)(z,n-s) = ((int)(z%(int)(pow(2,s)))) / (int)(pow(2,s-1));
    }
  }
}

vector<double> Original::getp()
{
  return _p1;
}

double Original::__getprop(vector<double>& p, int feat, int ind)
{
  // assert(ind < _sizeAlphabet);
  // assert(feat < _features.size());
  bool   found       = true;
  double sum         = 0.0;
  int    featuresize = 0;
  for(int i = 0; i < _sizeAlphabet; i++)
  {
    featuresize = _features[feat].size();
    for(int j = 0; j < featuresize; j++)
    {
      if((*_alphabet)(i,_features[feat][j]) != (*_alphabet)(ind,_features[feat][j]))
      {
        found = false;
        break;
      }
    }
    if(found)
    {
      sum += p[i];
    }
    found = true;
  }
  return sum;
}

double Original::getMarginalProp(int ind, vector<int> feat, vector<double> p)
{
  double sum = 0.0;
  bool found = true;
  for(int i = 0; i < _sizeAlphabet; i++)
  {
    for(int j = 0; j < feat.size(); j++)
    {
      if((*_alphabet)(i,feat[j]) != (*_alphabet)(ind,feat[j]))
      {
        found = false;
        break;
      }
    }
    if(found)
    {
      sum += p[i];
    }
    found = true;
  }
  return sum;
}
//p(featMarg|featCond)
double Original::getConditionalProp(vector<int> featMarg, vector<int> featCond,
                                    int ind, vector<double> p)
{
  double sum = 0.0;
  featMarg.insert(featMarg.end(), featCond.begin(), featCond.end());
  sum = getMarginalProp(ind, featMarg, p) / getMarginalProp(ind, featCond, p);
  return sum;
}

double Original::calculateKL(int iterations)
{
  double sum = 0.0;
  if(iterations == 0)
  {
#pragma omp parallel for
    for(int i = 0;i < _sizeAlphabet; i++)
    {
      if(_p1[i] > 0)
      {
        sum += _p1[i] * (log2(_p1[i]) - log2(_p2[i]));
      }
    }
  }
  else
  {
#pragma omp parallel for
    for(int i = 0; i < _sizeAlphabet; i++)
    {
      if(_p2[i] > 0)
      {
        sum += _p2[i] * (log2(_p2[i]) - log2(_p1[i]));
      }
    }
  }
  return sum;
}

double Original::calculateKL(vector<double> p, vector<double> q)
{
  double sum = 0.0;
  for(int i = 0; i < _sizeAlphabet; i++)
  {
    if(p[i] > 0)
    {
      sum += p[i] * (log2(p[i]) - log2(q[i]));
    }
  }
  return sum;
}

double Original::calculateConditionalKL(vector<double> p, vector<double> q,
                                        vector<int> featMarg, vector<int> featCond)
{
  double sum = 0.0;
  for(int i = 0; i < _sizeAlphabet; i++)
  {
    if(p[i] > 0)
    {
      sum += p[i] * (log2(getConditionalProp(featMarg,featCond,i,p))
                     - log2( getConditionalProp(featMarg,featCond,i,q)));
    }
  }
  return sum;
}

void Original::iterate(double threshold)
{
  int iterations = 0;
  double kl      = 2.0 * threshold;
  while( kl > threshold || iterations <= _features.size())
  {
    iterations++;
    cout << "Iteration: " << iterations << endl;
    int featIndex = (iterations-1) % _features.size();
    if((iterations % 2) != 0)
    {
#pragma omp parallel for
      for(int i = 0; i < _sizeAlphabet; i++)
      {
        _p2[i] = __getprop(_targetp, featIndex, i) * _p1[i];
        if(_p2[i] != 0)
        {
          _p2[i] = _p2[i] / __getprop(_p1,featIndex,i);
        }
      }
    }
    else
    {
#pragma omp parallel for
      for(int i = 0; i < _sizeAlphabet; i++)
      {
        _p1[i] = __getprop(_targetp, featIndex,i) * _p2[i];
        if(_p1[i] != 0)
        {
          _p1[i] = _p1[i] / __getprop(_p2,featIndex,i);
        }
      }
    }
    kl = calculateKL(iterations % 2);
    cout << "  KL: " << kl << endl;
  }
}

void Original::iterate(int iterations)
{
  for(int it = 1; it <= iterations; it++)
  {
    int featIndex = (it-1) % _features.size();
    if((it % 2) != 0)
    {
#pragma omp parallel for
      for(int i = 0; i < _sizeAlphabet; i++)
      {
        _p2[i] = __getprop(_targetp, featIndex,i) * _p1[i];
        if(_p2[i] != 0)
        {
          _p2[i] = _p2[i] / __getprop(_p1,featIndex,i);
        }
      }
    }
    else
    {
#pragma omp parallel for
      for(int i = 0; i < _sizeAlphabet; i++)
      {
        _p1[i] = __getprop(_targetp, featIndex,i) * _p2[i];
        if(_p1[i] != 0)
        {
          _p1[i] = _p1[i] / __getprop(_p2,featIndex,i);
        }
      }
    }
  }
}

