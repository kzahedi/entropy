#ifndef __SCGISGP_H__
#define __SCGISGP_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>
#include <entropy++/Matrix.h>
#include "Feature.h"
#include "InstanceMatrix.h"

using namespace std;

class SCGISgp{
public:
	SCGISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test);
	~SCGISgp();
	double 	prop(int Feati,int Featj,double ValX,double ValY);
	double 	prop(int rowX,vector<vector<double> > Y, int rowY);
	double 	prop(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY);
	double	getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY);
	double 	getconv(int i);
	int    	getsizeconv();

private:
	void			__scgis(int maxit, double konv,bool test, double lambdadeltaval,double sigma);
	double**** 		__getobs();

	int 			_sizeX;
	int 			_sizeY;
	int 			_sizeColValX;
	int 			_sizeColValY;
	int 			_sizeRowValX;
	int 			_sizeRowValY;
	double**** 		_observed;
	double****		_expected;
	double****		_exponent;
	double***		_normaliser;
	double****		_delta;
	vector<double> 	_conv;

	InstanceMatrix 	*_FM;
	DContainer 		*_Y;
	DContainer 		*_X;
	DContainer 		*_valY;
	DContainer 		*_valX;
};


#endif
