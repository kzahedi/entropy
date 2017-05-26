#include <entropy++/iterativescaling/GIS.h>
#include <entropy++/defs.h>

// #include <omp.h>

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
  createUniqueContainer();
  countObservedFeatures();
  _s.resize(Yalphabet->rows());

  _deltaMatcher = new DeltaMatcher(Xdata->rows());
  for(int x = 0; x < Xdata->rows(); x++)
  {
    vector<unsigned long> x_row = Xdata->row(x);
    // for(int y = 0; y < Yalphabet->rows(); y++)
    // for(int y = 0; y < Ydata->rows(); y++)
    {
      // vector<unsigned long> y_row = Yalphabet->row(y);
      vector<unsigned long> y_row = Ydata->row(x);
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
  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    (*d)->setExpected(0.0);
  }

#pragma omp parallel for
  for(int j = 0; j < Xdata->rows(); j++)
  {
    vector<unsigned long> x_row = Xdata->row(j);

    double *s = new double[Yalphabet->rows()];
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      // cout << omp_get_thread_num() << " " << omp_get_num_threads() << " " << y << endl;
      s[y] = 0.0;
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(vector<Delta*>::const_iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
          s[y] += (*d)->lambda();
        }
      }
    }

    double z = 0.0;
    for(int i = 0; i < Yalphabet->rows(); i++) z += exp(s[i]);

    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(vector<Delta*>::iterator d = _deltaMatcher->begin(j); d != _deltaMatcher->end(j); d++)
      {
        if((*d)->matchY(y_row))
        {
          (*d)->setExpected((*d)->expected() + exp(s[y]) / z);
        }
      }
    }
    delete[] s;
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
        if((*d)->matchXY(x_row, y_row)) f += 1.0;
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

  // cout << "**********" << endl;
  // for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  // {
    // cout << **d << endl;
  // }

}

double GIS::error()
{
  return _error;
}
