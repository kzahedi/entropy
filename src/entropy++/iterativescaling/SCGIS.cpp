#include <entropy++/iterativescaling/SCGIS.h>

#include <iostream>

using namespace std;

using namespace entropy::iterativescaling;

SCGIS::SCGIS() : Model()
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
  _error = 0.0;
  cout << deltas.size() << " * " << Yalphabet->rows() << " * " << Xdata->rows() << endl;
  cout << "Iterations: " << (deltas.size() * 2.0 * Yalphabet->rows() * Xdata->rows()) << endl;
  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    (*d)->setExpected(0);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(int j = 0; j < Xdata->rows(); j++)
      {
        vector<unsigned long> x_row = Xdata->row(j);
        if((*d)->matchXY(x_row, y_row))
        {
          (*d)->setExpected((*d)->expected() + exp((*_s)(j,y)) / (*_z)(j,0));
        }
      }
    }
    delta = log((*d)->observed() / (*d)->expected());
    _error += fabs((*d)->observed() - (*d)->expected());
    (*d)->setLambda((*d)->lambda() + delta);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(int j = 0; j < Xdata->rows(); j++)
      {
        vector<unsigned long> x_row = Xdata->row(j);
        if((*d)->matchXY(x_row, y_row))
        {
          (*_z)(j,0) -= exp((*_s)(j,y));
          (*_s)(j,y) += delta;
          (*_z)(j,0) += exp((*_s)(j,y));
        }
      }
    }
  }
}

double SCGIS::error()
{
  return _error;
}
