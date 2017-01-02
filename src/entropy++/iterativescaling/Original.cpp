#include "Original.h"
#include <entropy++/powi.h>

#include <iostream>

using namespace entropy::iterativescaling;
using namespace std;

Original::Original()
{
  _joint       = NULL;
  _conditional = NULL;
  _marginal    = NULL;
  _X           = NULL;
  _n           = 0;
  _initCalled  = false;
}

Original::~Original()
{
  if(_X != NULL) delete _X;
}

void Original::setData(ULContainer *X)
{
  if(_X->max() > 1)
  {
    cerr << "This implementation currently only works with binary data" << endl;
    exit(-1);
  }
  _X = X;
}

void Original::addFeature(vector<int> f)
{
  _features.push_back(f);
}

void Original::iterate()
{
  if(_initCalled == false)
  {
    cerr << "Init not called yet" << endl;
    exit(-1);
  }
  _n++;
}

void Original::init()
{
  _initCalled = true;
  _marginal   = new Matrix*[_features.size()];

  _Xalphabet = _X->unique();
  _emperical = new Matrix(1, powi(_Xalphabet->rows(), 2));
  for(int r = 0; r < _X->rows(); r++)
  {
    int index = 0;
    for(int c = 0; c < _X->columns(); c++)
    {
      index += _X->get(r,c) * powi(2,c);
    }
    (*_emperical)(0, index) += 1.0;
  }
  (*_emperical) /= (double)(_X->rows());

  for(int i = 0; i < (int)_features.size(); i++)
  {
    for(int r = 0; r < _X->rows(); r++)
    {
      int index = 0;
      for(vector<int>::iterator f = _features[i].begin(); f != _features[i].end(); f++)
      {

      }
    }
  }
}
