#ifndef _GIS_H_
#define _GIS_H_
#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"
#include "FeatureMatrix.h"

class GIS {
public:
	GIS(int ColValY,DContainer &eX,DContainer &aX, DContainer &aY);
	GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test,bool time,int seconds);
	~GIS();
	double 	prop(int rowX,vector<vector<double> > Y, int rowY);
	double 	prop(int Feati,int Featj,double ValX,double ValY);
	double 	propm(vector<vector<double> > X,int rowX,vector<vector<double> > Y);
	double 	propm(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY);
	void 	setFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY,double valuelambda);
	double	getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY);
	double 	getconv(int i);
	int    	getsizeconv();
	int		getIterations();

private:
	double**** 		__getobs();
	double   		__getFeatconst();
	void 			__getexp();
	void		 	__gis(int maxit, double konv, bool test);
	void 			__gis(int maxit, double konv, bool test,int seconds);

	double****		_expected;
	double**** 		_observed;
	double*** 		_exponent;
	double** 		_normaliser;
	vector<double> 	_conv;
	int 			_sizeX;
	int 			_sizeY;
	int 			_sizeColValX;
	int 			_sizeColValY;
	int 			_sizeRowValX;
	int 			_sizeRowValY;
	int				_iterations;

	FeatureMatrix 	*_FM;
	DContainer 		*_Y;
	DContainer 		*_X;
	DContainer 		*_valY;
	DContainer 		*_valX;
};

#endif
