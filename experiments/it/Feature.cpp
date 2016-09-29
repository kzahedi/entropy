#include "Feature.h"

Feature::Feature()
{
       _X = new DContainer(0, 0);
       _Y = new DContainer(0, 0);
   _sizeY = _Y->rows();
   _sizeX = _X->rows();
  _lambda = new Matrix(_sizeX,_sizeY);
  _sizeDeltaX = 0;
  _sizeDeltaY = 0;

}

//(Eingabealphabet, Anzahl der X und Y Werte, Anzahl der Testwerte, ein Startwert fuer alle  lambda)
Feature::Feature(DContainer &aX, DContainer &aY,int colValX, int colValY, int systXsize,int systYsize , double valuelambda)
{
  _X                    = &aX;
  _Y                    = &aY;
  assert(_X->columns() == 1);
  assert(_Y->columns() == 1);
  _sizeDeltaX           = pow(_X->rows(),systXsize);
  _sizeDeltaY           = pow(_Y->rows(),systYsize);
  _sizeY                = _Y->rows();
  _sizeX                = _X->rows();
  _lambda               = new Matrix(_sizeDeltaX,_sizeDeltaY);
  for(int i=0; i< _sizeDeltaX; i++)
  {
    for(int j=0; j< _sizeDeltaY; j++)
    {
      (*_lambda)(i,j) = valuelambda;
    }
  }
}
//alle lambda explizit ueber die Matrix setzen
Feature::Feature(DContainer &aX, DContainer &aY,int colValX, int colValY, int systXsize,int systYsize, Matrix &lambda)
{
  assert(aX.columns()==1);
  assert(aY.columns()==1);

  _sizeDeltaX  = pow(aX.rows(),systXsize);
  _sizeDeltaY  = pow(aY.rows(),systYsize);
  assert(lambda.rows()==_sizeDeltaX);
  assert(lambda.cols()==_sizeDeltaY);
  _lambda = &lambda;
  _X      = &aX;
  _Y      = &aY;
  _sizeY  = _Y->rows();
  _sizeX  = _X->rows();
}

Feature::~Feature()
{
  delete _lambda;
  delete _X;
  delete _Y;
}

double Feature::getlambda(int i, int j)
{
  assert(i < _sizeDeltaX && j < _sizeDeltaY);
  return (*_lambda)(i,j);
}

void Feature::setlambda(int i, int j, double newvalue)
{
  assert(i < _sizeDeltaX && j < _sizeDeltaY);
  (*_lambda)(i,j) = newvalue;
}

Feature& Feature::operator=(const Feature& c)
{
  this-> _sizeX = c._sizeX;
  this-> _sizeY = c._sizeY;
  this-> _sizeDeltaX = c._sizeDeltaX;
  this-> _sizeDeltaY = c._sizeDeltaY;
  this-> _X     = new DContainer(_sizeX,1);
  this-> _Y     = new DContainer(_sizeY,1);

  _Y            = c._Y;
  _X            = c._X;

  _lambda       = new Matrix(_sizeX,_sizeY);
  this->_lambda = c._lambda;

  return *this;
}

