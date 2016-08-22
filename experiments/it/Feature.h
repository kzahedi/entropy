#ifndef __FEATURE_H__
#define __FEATURE_H__

#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

class Feature
{
public:
	//Konstruktoren:
	 Feature();

	Feature(DContainer &aX, DContainer &aY, double valuelambda){ //(Eingabealphabet, Startwert fuer lambda)
		X= &aX;
		Y= &aY;
		assert((*X).columns()==1);
		assert((*Y).columns()==1);
		sizeY= (*Y).rows();
		_sizeX= (*X).rows();
		lambda= new double*[_sizeX];
		for(int m=0; m<_sizeX; m++)
			  lambda[m]= new double[sizeY];

		for(int i=0; i< _sizeX; i++){
			for(int j=0; j< sizeY; j++){
				lambda[i][j]=valuelambda;
			}
		}
	}
	Feature(bool binaer,double valuelambda){ //Eingabealphabete nur 1, -1
		X= new DContainer(2,1);
		*X << 1 << -1;
		Y= new DContainer(2,1);
		*Y << 1 << -1;
		sizeY= (*Y).rows();
		_sizeX= (*X).rows();
		lambda= new double*[_sizeX];
				for(int m=0; m<_sizeX; m++)
					  lambda[m]= new double[sizeY];

				for(int i=0; i< _sizeX; i++){
					for(int j=0; j< sizeY; j++){
						lambda[i][j]=valuelambda;
					}
				}
	}
	Feature(DContainer &aX, DContainer &aY,  double** lambda1){
		lambda= lambda1;
		X= &aX;
		Y= &aY;
		sizeY= (*Y).rows();
		_sizeX= (*X).rows();
		assert((*X).columns()==1);
		assert((*Y).columns()==1);
		//Groesse von lambda1 testen
		assert(lambda[_sizeX-1][sizeY-1]&& !lambda[_sizeX][sizeY]);  //
	}
	~Feature(){

	}

	double getlambda(int i, int j);
	void setlambda(int i, int j, double newvalue );

	double value(double x,double y){
		double val=0;
		for(int i=0; i< _sizeX; i++){
			for(int j=0; j< sizeY; j++){
				double a=(*X).get(i,0);
				double b=(*Y).get(j,0);
				val= val + lambda[i][j]*delta(a, b, x, y);
			}
		}
		return val;
	}

	Feature& operator=(const Feature& c){
		this-> _sizeX= c._sizeX;
		this-> sizeY= c.sizeY;
		this-> X = new DContainer(_sizeX,1);
		this-> Y = new DContainer(sizeY,1);

		Y=c.Y;
		X=c.X;

		this-> lambda = new double*[_sizeX];
		for(int m=0; m<_sizeX; m++)
			  lambda[m]= new double[sizeY];
		for(int i=0; i< _sizeX; i++){
						for(int j=0; j< sizeY; j++){
							this->lambda[i][j]=c.lambda[i][j];
						}
		}
		return *this;
	}


private:

	int delta(double ax, double ay, double x, double y){     // alphabet x,y und eingegebenes x,y
		if(ax== x && ay== y){
			return 1;
		}
		else return 0;
	}


	int _sizeX;
	int sizeY;
	DContainer *X;
	DContainer *Y;
	double** lambda;


};

#endif
