#include <entropy++/iterativescaling/GIS.h>

using namespace entropy::iterativescaling;

GIS::GIS() : Model()
{
}

void GIS::init()
{
  createUniqueContainer();
  countObservedFeatures();
}

void GIS::iterate()
{
  generateExpected();

  _error = 0.0;

  double max = 1.0;

  // get f^#
  for(int j = 0; j < Xdata->rows(); j++)
  {
    vector<unsigned long> x_row = Xdata->row(j);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> Y = Yalphabet->row(y);
      double f = 0.0;
      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, Y)) f += 1.0;
      }
      if(f > max) max = f;
    }
  }

  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    double n = (*d)->lambda() + fabs(1.0/max) * log((*d)->observed() / (*d)->expected());
    (*d)->setLambda(n);
    _error += fabs((*d)->observed() - (*d)->expected());
  }
}

double GIS::error()
{
  return _error;
}
