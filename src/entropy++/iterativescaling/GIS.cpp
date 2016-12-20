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

  int index = 0;

  _error = 0.0;

  double max = 0.0;
  double n   = 0.0;
  double o   = 0.0;


  // get f^#
  for(int j = 0; j < Xdata->rows(); j++)
  {
    vector<unsigned long> x_row = Xdata->row(j);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> Y = Yalphabet->row(y);

      double f = 0.9;

      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, Y))
        {
          f += (*d)->lambda();
        }
      }
      if(f > max) max = f;
  }

    for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
    {
      // double o = (*d)->lambda(); // old
      double n = (*d)->lambda() + (1.0/max) * log((*d)->observed() / (*d)->expected());
      (*d)->setLambda(n);
      _error += fabs((*d)->observed() - (*d)->expected());
    }
  }
}

double GIS::error()
{
  return _error;
}
