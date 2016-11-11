#include "Feature.h"

using namespace entropy::iterativescaling;

Feature::Feature()
{
       _X = new DContainer(0, 0);
       _Y = new DContainer(0, 0);
   _sizeY = _Y->rows();
   _sizeX = _X->rows();
  _lambda = new SparseMatrix();
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
  _lambda               = new SparseMatrix(valuelambda);
}
//alle lambda explizit ueber die Matrix setzen
Feature::Feature(DContainer &aX, DContainer &aY,int colValX, int colValY, int systXsize,int systYsize, SparseMatrix &lambda)
{
  assert(aX.columns()==1);
  assert(aY.columns()==1);

  _sizeDeltaX  = pow(aX.rows(),systXsize);
  _sizeDeltaY  = pow(aY.rows(),systYsize);
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
int Feature::getLambdaSize(){
	int j= _lambda->size();
	return j;
}
double Feature::getLambda(int i, int j)
{
  assert(i < _sizeDeltaX && j < _sizeDeltaY);
  return (*_lambda)(i,j);
}

void Feature::setLambda(int i, int j, double newvalue)
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

  _lambda       = new SparseMatrix();
  this->_lambda = c._lambda;

  return *this;
}

