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
  _systX = NULL;
  _systY = NULL;

}

//(Eingabealphabet, Startwert fuer lambda)
Feature::Feature(DContainer &aX, DContainer &aY,vector<int>& systX,vector<int>& systY, double valuelambda)
{
	cout << "hier 1 " << endl;
  _X                    = &aX;
  _Y                    = &aY;
  assert(_X->columns() == 1);
  assert(_Y->columns() == 1);
  _systX                = &systX;
  _systY                = &systY;
  _sizeDeltaX           = pow(_X->rows(),_systX->size());
  _sizeDeltaY           = pow(_Y->rows(),_systY->size());
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

Feature::Feature(DContainer &aX, DContainer &aY,vector<int>& systX,vector<int>& systY, Matrix &lambda)
{
  assert(aX.columns()==1);
  assert(aY.columns()==1);
  _systX = &systX;
  _systY = &systY;
  _sizeDeltaX  = pow(aX.rows(),_systX->size());
  _sizeDeltaY  = pow(aY.rows(),_systY->size());
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

double Feature::value(vector<double> x,vector<double> y)
{
  assert(x.size()==_systX->size() && y.size()== _systY->size());
  double val=0.0;
  for(int i=0; i< _sizeDeltaX; i++)
  {
    for(int j=0; j< _sizeDeltaY; j++)
    {
      val+= (*_lambda)(i,j)*delta(i, j, x, y);
    }
  }
  return val;
}
vector<double> Feature::index(int index,bool x)
{
  int sizeAlph;
  int sizeSyst;
  vector<double> zeile;
  double z;
  if(x){
	  sizeAlph=_sizeX;
	  sizeSyst=_systX->size();
  }
  else{
	  sizeAlph=_sizeY;
	  sizeSyst=_systY->size();
  }
  for(int i=sizeSyst; i>0 ;i--)
  {
	  z = index % (int) (pow(sizeAlph,i));
	  z = z/ (pow(sizeAlph,i-1));
	  int j= (int) z;
	  if(x){
		  zeile.push_back((*_X)(j,0));
	  }
	  else{
		  zeile.push_back((*_Y)(j,0));
	  }
  }
  return zeile;
}
Feature& Feature::operator=(const Feature& c)
{
  this-> _sizeX = c._sizeX;
  this-> _sizeY = c._sizeY;
  this-> _X     = new DContainer(_sizeX,1);
  this-> _Y     = new DContainer(_sizeY,1);

  _Y            = c._Y;
  _X            = c._X;

  _lambda       = new Matrix(_sizeX,_sizeY);
  this->_lambda = c._lambda;

  return *this;
}
// alphabet ax,ay und eingegebenes x,y
int Feature::delta(int indexX, int indexY, vector<double> x, vector<double> y)
{
  vector<double> aX = index(indexX,true);
  vector<double> aY = index(indexY,false);
  assert( aX.size() == x.size() && aY.size() == y.size());
  bool equ = true;
  for(int i=0; i<_systX->size();i++)
  {
	  if(aX[i] != x[i]){
		  equ=false;
	  }
  }
  for(int j=0; j<_systY->size(); j++)
  {
	  if(aY[j] != y[j]){
		  equ=false;
	  }
  }
  if(equ){
	  return 1;
  }
  else{
	  return -1;
  }
}
