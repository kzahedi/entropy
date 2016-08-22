#include "Feature.h"

Feature::Feature() {
	_X = new DContainer(0, 0);
	_Y = new DContainer(0, 0);
	_sizeY = (*_Y).rows();
	_sizeX = (*_X).rows();
	_lambda = new double*[_sizeX];
	for (int m = 0; m < _sizeX; m++)
		_lambda[m] = new double[_sizeY];
}

//(Eingabealphabet, Startwert fuer lambda)
Feature::Feature(DContainer &aX, DContainer &aY, double valuelambda){
	_X= &aX;
	_Y= &aY;
	assert((*_X).columns()==1);
	assert((*_Y).columns()==1);
	_sizeY= (*_Y).rows();
	_sizeX= (*_X).rows();
	_lambda= new double*[_sizeX];
	for(int m=0; m<_sizeX; m++)
		  _lambda[m]= new double[_sizeY];

	for(int i=0; i< _sizeX; i++){
		for(int j=0; j< _sizeY; j++){
			_lambda[i][j]=valuelambda;
		}
	}
}
//Eingabealphabete nur 1, -1
Feature::Feature(bool binaer,double valuelambda){
		_X= new DContainer(2,1);
		*_X << 1 << -1;
		_Y= new DContainer(2,1);
		*_Y << 1 << -1;
		_sizeY= (*_Y).rows();
		_sizeX= (*_X).rows();
		_lambda= new double*[_sizeX];
				for(int m=0; m<_sizeX; m++)
					  _lambda[m]= new double[_sizeY];

				for(int i=0; i< _sizeX; i++){
					for(int j=0; j< _sizeY; j++){
						_lambda[i][j]=valuelambda;
					}
				}
	}
Feature::Feature(DContainer &aX, DContainer &aY,  double** lambda){
	_lambda= lambda;
	_X= &aX;
	_Y= &aY;
	_sizeY= (*_Y).rows();
	_sizeX= (*_X).rows();
	assert((*_X).columns()==1);
	assert((*_Y).columns()==1);
}

double Feature::getlambda(int i, int j) {
	assert(i < _sizeX && j < _sizeY);
	return _lambda[i][j];
}

void Feature::setlambda(int i, int j, double newvalue)
{
	assert(i < _sizeX && j < _sizeY);
	_lambda[i][j] = newvalue;
}
double Feature::value(double x,double y){
		double val=0;
		for(int i=0; i< _sizeX; i++){
			for(int j=0; j< _sizeY; j++){
				double a=(*_X).get(i,0);
				double b=(*_Y).get(j,0);
				val= val + _lambda[i][j]*__delta(a, b, x, y);
			}
		}
		return val;
	}
Feature& Feature::operator=(const Feature& c){
	this-> _sizeX= c._sizeX;
	this-> _sizeY= c._sizeY;
	this-> _X = new DContainer(_sizeX,1);
	this-> _Y = new DContainer(_sizeY,1);

	_Y=c._Y;
	_X=c._X;

	this->_lambda = new double*[_sizeX];
	for(int m=0; m<_sizeX; m++)
		  this->_lambda[m]= new double[_sizeY];
	for(int i=0; i< _sizeX; i++){
					for(int j=0; j< _sizeY; j++){
						this->_lambda[i][j]=c._lambda[i][j];
					}
	}
	return *this;
}
int Feature::__delta(double ax, double ay, double x, double y){     // alphabet x,y und eingegebenes x,y
	if(ax== x && ay== y){
		return 1;
	}
	else return 0;
}
//  #include <tuple>
// typedef std::tuple<int,int> Tuple;
// Tuple** matrix;

