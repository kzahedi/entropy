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

class GISgp : public IT{

public:
	GISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,double lambdadeltaval, double sigma,int maxit,double konv, bool test,bool time,int seconds);
	~GISgp();
	double 	getconv(int i);
	int    	getsizeconv();
	int		getIterations();

private:
	double   		__getFeatconst();
	void 			__getexp();
	void 			__gisgp(int maxit, double konv, double lambdadelta, double sigma,bool test);
	void 			__gisgp(int maxit, double konv, double lambdadelta, double sigma,bool test,int seconds);

	double****		_expected;
	double****		_delta;
	double* 		_exponent;
	double			_normaliser;
	vector<double> 	_conv;
	int 			_iterations;

};
#endif

