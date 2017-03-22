#include <entropy++/iterativescaling/SCGIS.h>
#include <iostream>
#include <omp.h>

// #include <glog/logging.h>

using namespace std;

using namespace entropy::iterativescaling;

SCGIS::SCGIS() : IS()
{
  _s = NULL;
  _z = NULL;
}

SCGIS::~SCGIS()
{
  if(_z != NULL) delete _z;
  if(_s != NULL) delete _s;
}

void SCGIS::init()
{
  createUniqueContainer();
  countObservedFeatures();
  _z = new Matrix(Xdata->rows(), 1, Yalphabet->rows());
  _s = new Matrix(Xdata->rows(), Yalphabet->rows(), 0.0);
}

void SCGIS::iterate()
{
  double delta = 0.0;
  double e     = 0.0;
  _error       = 0.0;

#pragma omp parallel for private(delta),shared(deltas,e)
  for(int i = 0; i < (int)deltas.size(); i++)
  {
    Delta *d = deltas[i];
    d->setExpected(0);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(int j = 0; j < Xdata->rows(); j++)
      {
        vector<unsigned long> x_row = Xdata->row(j);
        if(d->matchXY(x_row, y_row))
        {
          d->setExpected(d->expected() + exp((*_s)(j,y)) / (*_z)(j,0));
        }
      }
    }
    delta   = log(d->observed() / d->expected());
    e       = d->observed() - d->expected();
    _error += e * e;
    d->setLambda(d->lambda() + delta);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(int j = 0; j < Xdata->rows(); j++)
      {
        vector<unsigned long> x_row = Xdata->row(j);
        if(d->matchXY(x_row, y_row))
        {
          (*_z)(j,0) -= exp((*_s)(j,y));
          (*_s)(j,y) += delta;
          (*_z)(j,0) += exp((*_s)(j,y));
        }
      }
    }
  }
  _error = sqrt(_error);
}

double SCGIS::error()
{
  return _error;
}
