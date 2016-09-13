#ifndef __GISGP_H__
#define __GISGP_H__

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
#include "SCGIS.h"
#include "GIS.h"

using namespace std;

class GISgp{

public:
	GISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,double lambdadeltaval, double sigma,int maxit,double konv, bool test);
	~GISgp();
	double 	prop(int rowX,vector<vector<double> > Y, int rowY);
	double 	prop(int Feati,int Featj,double ValX,double ValY);
	double 	prop(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY);
	void 	setFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY,double valuelambda);
	double	getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY);
	double 	getconv(int i);
	int    	getsizeconv();

private:
	double**** 		__getobs();
	double   		__getFeatconst();
	void 			__getexp();
	void 			__gisgp(int maxit, double konv, double lambdadelta, double sigma,bool test);

	double****		_expected;
	double**** 		_observed;
	double****		_lambdadelta;
	double*** 		_exponent;
	double** 		_normaliser;
	vector<double> 	_conv;
	int 			_sizeX;
	int 			_sizeY;
	int 			_sizeColValX;
	int 			_sizeColValY;
	int 			_sizeRowValX;
	int 			_sizeRowValY;

	FeatureMatrix 	*_FM;
	DContainer 		*_Y;
	DContainer 		*_X;
	DContainer 		*_valY;
	DContainer 		*_valX;
};
#endif
