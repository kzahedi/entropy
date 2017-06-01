#include <entropy++/iterativescaling/SCGIS.h>
#include <entropy++/defs.h>
#include <iostream>
#ifdef USE_OPENMP
#  include <omp.h>
#endif // USE_OPENMP
#include <math.h>

#define EPSILON 0.00001

#include <glog/logging.h>

using namespace std;

using namespace entropy::iterativescaling;

SCGIS::SCGIS() : IS()
{
  _s          = NULL;
  _z          = NULL;
  _rowMatcher = NULL;
}

SCGIS::~SCGIS()
{
  if(_z          != NULL) delete _z;
  if(_s          != NULL) delete _s;
  if(_rowMatcher != NULL) delete _rowMatcher;
}

void SCGIS::init()
{
  createUniqueContainer();
  countObservedFeatures();
  _z          = new Matrix(Xdata->rows(), 1, Yalphabet->rows());
  _s          = new Matrix(Xdata->rows(), Yalphabet->rows(), 0.0);
  _rowMatcher = new RowMatcher(deltas.size());

#pragma omp parallel for
  for(int i = 0; i < (int)deltas.size(); i++)
  {
    Delta *d = deltas[i];
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      if(d->matchY(y_row)) _rowMatcher->add_y(i, y);
    }
    for(int x = 0; x < Xdata->rows(); x++)
    {
      vector<unsigned long> x_row = Xdata->row(x);
      if(d->matchX(x_row)) _rowMatcher->add_x(i, x);
    }
  }
}

void SCGIS::iterate()
{
  _error = 0.0;

  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    (*d)->setExpected(0);
  }

// #pragma omp parallel for shared(deltas,e,_rowMatcher)
  for(int i = 0; i < (int)deltas.size(); i++)
  {
    Delta *d = deltas[i];
    // for(vector<int>::iterator y = _rowMatcher->y_begin(i); y != _rowMatcher->y_end(i); y++)
#pragma omp parallel for
    for(int yi = 0; yi < _rowMatcher->y_size(i); yi++)
    {
      int y = _rowMatcher->y(i, yi);
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(vector<int>::iterator j = _rowMatcher->x_begin(i); j != _rowMatcher->x_end(i); j++)
      {
        vector<unsigned long> x_row = Xdata->row(*j);
        if(d->matchXY(x_row, y_row))
        {
          d->updateExpected(exp((*_s)(*j,y)) / (*_z)(*j,0));
          if(isnan(d->expected()) || isinf(d->expected()))
          {
            cerr << (*_s)(*j,y) << " " << exp((*_s)(*j,y)) << "/" << (*_z)(*j,0) << " " <<  exp((*_s)(*j,y)) / (*_z)(*j,0) << endl;
          }
          // must do better
          if(d->expected() < EPSILON) d->setExpected(EPSILON);
        }
      }
    }
    double delta = log(d->observed() / d->expected());
    double e     = d->observed() - d->expected();

    _error += e * e;

    if(isnan(_error))
    {
      cerr << "delta: " << delta << " = log(" << d->observed() << " / " << d->expected() << ")" << endl;
      exit(-1);
    }
    d->updateLambda(delta);
    // for(vector<int>::iterator y = _rowMatcher->y_begin(i); y != _rowMatcher->y_end(i); y++)
    // {
#pragma omp parallel for
    for(int yi = 0; yi < _rowMatcher->y_size(i); yi++)
    {
      int y = _rowMatcher->y(i, yi);
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(vector<int>::iterator j = _rowMatcher->x_begin(i); j != _rowMatcher->x_end(i); j++)
      {
        vector<unsigned long> x_row = Xdata->row(*j);
        if(d->matchXY(x_row, y_row))
        {
#pragma omp critical
          {
            (*_z)(*j,0)  -= exp((*_s)(*j,y));
            (*_s)(*j,y) += delta;
            (*_z)(*j,0)  += exp((*_s)(*j,y));
          }
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
