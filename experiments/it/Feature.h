#ifndef __FEATURE_H__
#define __FEATURE_H__

#include <entropy++/Container.h>

#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>
#include <vector>

using namespace std;

class Feature
{
public:

	Feature();
	Feature(DContainer &aX, DContainer &aY, double valuelambda);
	Feature(DContainer &aX, DContainer &aY,  double** lambda);
	Feature(bool binaer,double valuelambda);
	~Feature();

	friend std::ostream& operator<<(std::ostream& str,Feature& feature){
		str<< "Feature:" <<endl;
		for(int i=0; i<feature._sizeX; i++){
			for(int j=0; j<feature._sizeY; j++){
				str<< feature.getlambda(i,j) << " ";
			}
			str<< endl;
		}
	   return str;
	};

	double getlambda(int i, int j);
	void setlambda(int i, int j, double newvalue );
	double value(double x,double y);
	Feature& operator=(const Feature& c);


private:

	int __delta(double ax, double ay, double x, double y);

	int _sizeX;
	int _sizeY;
	DContainer *_X;
	DContainer *_Y;
	double** _lambda;

};

#endif
