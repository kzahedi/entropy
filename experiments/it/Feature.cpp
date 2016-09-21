#include "Feature.h"

Feature::Feature() {
	_X = new DContainer(0, 0);
	_Y = new DContainer(0, 0);
	_sizeY = _Y->rows();
	_sizeX = _X->rows();
	_lambda = new Matrix(_sizeX,_sizeY);
}

//(Eingabealphabet, Startwert fuer lambda)
Feature::Feature(DContainer &aX, DContainer &aY, double valuelambda){
	_X= &aX;
	_Y= &aY;
	assert(_X->columns()==1);
	assert(_Y->columns()==1);
	_sizeY=_Y->rows();
	_sizeX=_X->rows();
	_lambda = new Matrix(_sizeX,_sizeY);
	for(int i=0; i< _sizeX; i++){
		for(int j=0; j< _sizeY; j++){
			(*_lambda)(i,j)=valuelambda;
		}
	}
}
Feature::Feature(DContainer &aX, DContainer &aY, Matrix lambda){
	_lambda=&lambda;
	_X= &aX;
	_Y= &aY;
	_sizeY= _Y->rows();
	_sizeX= _X->rows();
	assert(_X->columns()==1);
	assert(_Y->columns()==1);
	assert(_lambda->rows()==(_sizeX));
	assert(_lambda->cols()==(_sizeY));
}
Feature:: ~Feature(){
	delete _lambda;
	delete _X;
	delete _Y;
}

double Feature::getlambda(int i, int j) {
	assert(i < _sizeX && j < _sizeY);
	return (*_lambda)(i,j);
}

void Feature::setlambda(int i, int j, double newvalue){
	assert(i < _sizeX && j < _sizeY);
	(*_lambda)(i,j) = newvalue;
}

double Feature::value(double x,double y){
		double val=0;
		for(int i=0; i< _sizeX; i++){
			for(int j=0; j< _sizeY; j++){
				double a=_X->get(i,0);
				double b=_Y->get(j,0);
				val+= (*_lambda)(i,j)*delta(a, b, x, y);
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

	_lambda = new Matrix(_sizeX,_sizeY);
	this->_lambda=c._lambda;

	return *this;
}
// alphabet ax,ay und eingegebenes x,y
int Feature::delta(double ax, double ay, double x, double y){
	if(ax== x && ay== y){
		return 1;
	}
	else return -1;
}


