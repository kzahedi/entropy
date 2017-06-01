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

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
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
  double e     = 0.0;
  _error       = 0.0;

#ifdef USE_OPENMP
#pragma omp parallel for shared(deltas,e,_rowMatcher)
#endif
  for(int i = 0; i < (int)deltas.size(); i++)
  {
#ifdef USE_OPENMP
    VLOG(100) << "working i: " << i << " thread: " << omp_get_thread_num() << " of " << omp_get_num_threads();
#endif
    Delta *d = deltas[i];
    d->setExpected(0);
    // for(int y = 0; y < Yalphabet->rows(); y++)
    for(vector<int>::iterator y = _rowMatcher->y_begin(i); y != _rowMatcher->y_end(i); y++)
    {
#ifdef USE_OPENMP
      VLOG(100) << "y: " << *y;
#endif
      vector<unsigned long> y_row = Yalphabet->row(*y);
      // for(int j = 0; j < Xdata->rows(); j++)
      for(vector<int>::iterator j = _rowMatcher->x_begin(i); j != _rowMatcher->x_end(i); j++)
      {
#ifdef USE_OPENMP
        VLOG(100) << "j: " << *j;
#endif
        vector<unsigned long> x_row = Xdata->row(*j);
        if(d->matchXY(x_row, y_row))
        {
          d->setExpected(d->expected() + exp((*_s)(*j,*y)) / (*_z)(*j,0));
          // must do better
          if(d->expected() < EPSILON) d->setExpected(EPSILON);
        }
      }
    }
    double delta   = log(d->observed() / d->expected());
    e       = d->observed() - d->expected();
    _error += e * e;
    if(isnan(_error))
    {
      cout << "delta: " << delta << " = log(" << d->observed() << " / " << d->expected() << ")" << endl;
    }
    d->setLambda(d->lambda() + delta);
    // for(int y = 0; y < Yalphabet->rows(); y++)
    for(vector<int>::iterator y = _rowMatcher->y_begin(i); y != _rowMatcher->y_end(i); y++)
    {
#ifdef USE_OPENMP
      VLOG(100) << "2 y: " << *y;
#endif
      vector<unsigned long> y_row = Yalphabet->row(*y);
      // for(int j = 0; j < Xdata->rows(); j++)
      for(vector<int>::iterator j = _rowMatcher->x_begin(i); j != _rowMatcher->x_end(i); j++)
      {
#ifdef USE_OPENMP
      VLOG(100) << "2 j: " << *j;
#endif
        vector<unsigned long> x_row = Xdata->row(*j);
        if(d->matchXY(x_row, y_row))
        {
          (*_z)(*j,0)  -= exp((*_s)(*j,*y));
          (*_s)(*j,*y) += delta;
          (*_z)(*j,0)  += exp((*_s)(*j,*y));
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
