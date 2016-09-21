#include "ITMatrix.h"

ITMatrix::ITMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda)
{
  _valX        = &eX;
  _valY        = &eY;
  _X           = &aX;
  _Y           = &aY;
  _sizeColValY = _valY->columns();
  _sizeColValX = _valX->columns();
  _sizeRowValX = _valX->rows();
  _sizeRowValY = _valY->rows();
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  _FA          = FeatureArray(valuelambda);
}

ITMatrix::ITMatrix()
{
  _valX        = new DContainer(0,0);
  _valY        = new DContainer(0,0);
  _X           = new DContainer(0,0);
  _Y           = new DContainer(0,0);
  _sizeColValY = 0;
  _sizeColValX = 0;
  _sizeRowValX = 0;
  _sizeRowValY = 0;
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  _FA          = FeatureArray(0);
}

ITMatrix::~ITMatrix()
{
  for(int i=0; i<_sizeColValX;i++)
  {
    delete _FA[i];
  }
  delete _FA;
}

double  ITMatrix::getFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY)
{
  assert(i<_sizeColValX && j<_sizeColValY);
  double lambda=_FA[i][j].getlambda(ilambdaX,ilambdaY);
  return lambda;
}

double ITMatrix::getFeatureArrayvalue(int i, int j,double ValX, double ValY)
{
  assert(i<_sizeColValX && j<_sizeColValY);
  double value=_FA[i][j].value(ValX,ValY);
  return value;
}

void  ITMatrix::setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda)
{
  assert(i<_sizeColValX && j<_sizeColValY);
  _FA[i][j].setlambda(ilambdaX,ilambdaY,valuelambda);
}

int   ITMatrix::getFeatureArraydelta(int i, int j,int idelta, int jdelta, double ValX, double ValY)
{
  assert(i<_sizeColValX && j<_sizeColValY);
  assert(idelta < _sizeX && jdelta < _sizeY);
  int delta=_FA[i][j].delta(_X->get(idelta,0),_Y->get(jdelta,0), ValX,ValY);
  return delta;
}

Feature** ITMatrix::FeatureArray( double valuelambda)
{
  Feature **FA;
  FA = new Feature*[_sizeColValX];
  for(int m=0; m< _sizeColValX; m++)
  {
    FA[m]= new Feature[_sizeColValY];
    for(int k=0;k<_sizeColValY;k++)
    {
      Feature *K=new Feature(*_X,*_Y,valuelambda);
      FA[m][k]=*K;
    }
  }
  return FA;
}
