#include <entropy++/iterativescaling/GIS.h>
#include <entropy++/defs.h>

#include <glog/logging.h>

#ifdef USE_OPENMP
#  include <omp.h>
#endif

// #include <glog/logging.h>

#define MIN_S 0.0000001

using namespace entropy::iterativescaling;

GIS::GIS() : IS()
{ }

GIS::~GIS()
{
  delete _deltaMatcher;
}

void GIS::init()
{
  VLOG(10) << "creating unique containers";
  createUniqueContainer();
  VLOG(10) << "counting observed features";
  countObservedFeatures();
  _s.resize(Yalphabet->rows());

  VLOG(10) << "Initialising DeltaMatcher";
  _deltaMatcher = new DeltaMatcher(Xdata->rows());
  for(int x = 0; x < Xdata->rows(); x++)
  {
    vector<unsigned long> x_row = Xdata->row(x);

    if (VLOG_IS_ON(10))
    {
      stringstream sst;
      sst << "x-row " << x << " data:";
      for(vector<unsigned long>::iterator i = x_row.begin(); i != x_row.end(); i++)
      {
        sst << " " << *i;
      }
      VLOG(10) << sst.str();
    }

    for(int y = 0; y < Yalphabet->rows(); y++)
    // for(int y = 0; y < Ydata->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      // vector<unsigned long> y_row = Ydata->row(x);

      if (VLOG_IS_ON(10))
      {
        stringstream sst;
        sst << "  y-row " << x << " data:";
        for(vector<unsigned long>::iterator i = y_row.begin(); i != y_row.end(); i++)
        {
          sst << " " << *i;
        }
        VLOG(10) << sst.str();
      }

      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
          bool found = false;
          for(vector<Delta*>::iterator dm = _deltaMatcher->begin(x); dm != _deltaMatcher->end(x); dm++)
          {
            if(*dm == *d)
            {
              found = true;
              break;
            }
          }
          if(found == false) _deltaMatcher->add(x, *d);
        }
      }
    }
  }
}

void GIS::iterate()
{
  VLOG(100) << "Setting expected to 0";
  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    (*d)->setExpected(0.0);
  }


#pragma omp parallel for
  for(int j = 0; j < Xdata->rows(); j++)
  {
    VLOG(10) << "Working on row " << j;
    vector<unsigned long> x_row = Xdata->row(j);

    vector<double> s;
    s.resize(Yalphabet->rows());
// #pragma omp parallel for
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      // cout << omp_get_thread_num() << " " << omp_get_num_threads() << " " << y << endl;
      s[y] = 0.0;
      vector<unsigned long> y_row = Yalphabet->row(y);
      // for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      for(vector<Delta*>::iterator d = _deltaMatcher->begin(j); d != _deltaMatcher->end(j); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
#pragma omp atomic update
          s[y] += (*d)->lambda();
        }
      }
    }

    double z = 0.0;
    for(int i = 0; i < Yalphabet->rows(); i++) z += exp(s[i]);

// #pragma omp parallel for
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      // ROW_OUTPUT("update x:", x_row);
      // ROW_OUTPUT("update y:", y_row);
      for(vector<Delta*>::iterator d = _deltaMatcher->begin(j); d != _deltaMatcher->end(j); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
          VLOG(100) << "found:  " << **d;
          VLOG(100) << "exp:    " << (*d)->expected();
          VLOG(100) << "z:      " << z;
          VLOG(100) << "exp(s): " << exp(s[y]);
          VLOG(100) << "delta:  " << exp(s[y])/z;
          // (*d)->setExpected((*d)->expected() + exp(s[y]) / z);
          (*d)->updateExpected(exp(s[y]) / z);
          VLOG(100) << "set to: " << **d;
        }
      }
    }
  } // j

  _error = 0.0;

  double max = 1.0;
  double e = 0.0;

  // get f^#
  for(int j = 0; j < Xdata->rows(); j++)
  {
    vector<unsigned long> x_row = Xdata->row(j);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      double f = 0.0;
      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
#pragma omp atomic update
          f += 1.0;
        }
      }
      if(f > max) max = f;
    }
  }

  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    double n = (*d)->lambda() + fabs(1.0/max) * log((*d)->observed() / (*d)->expected());
    (*d)->setLambda(n);
    e = (*d)->observed() - (*d)->expected();
    _error += e * e;
  }

  _error = sqrt(_error);

}

double GIS::error()
{
  return _error;
}
