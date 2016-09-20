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
#include "IT.h"

using namespace std;

class SCGISgp : public IT{
public:
	SCGISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,double lambdadeltaval, double sigma ,int maxit, double konv, bool test,bool time,int seconds);
	~SCGISgp();
	double 	getconv(int i);
	int    	getsizeconv();
	int		getIterations();

private:
	void			__scgis(int maxit, double konv,bool test, double lambdadeltaval,double sigma);
	void			__scgis(int maxit, double konv,bool test, double lambdadeltaval,double sigma,int seconds);

	int				_iterations;
	double****		_expected;
	double****		_exponent;
	double***		_normaliser;
	double****		_delta;
	vector<double> 	_conv;
};


#endif

