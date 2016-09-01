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
	GIS(int ColValY,DContainer &eX);
	GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double la, int maxit, double konv);
	GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test);
	double gis(int rowX,vector<vector<double> > Y, int rowY);
	double gis(int Feati,int Featj,double ValX,double ValY);
	void setFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY,double valuelambda);
	double getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY);
	double getconv(int i);
	int    getsizeconv();

private:
	double**** __getobs();
	void __gis(int maxit, double konv);
	double __getFeatconst();
	void __getexp(double**** &expect, double*** &exponent,double** &normaliser);
	vector<double> __gis(int maxit, double konv, bool test);

	vector<double> conv;
	int _sizeX;
	int _sizeY;
	int _sizeColValX;
	int _sizeColValY;
	int _sizeRowValX;
	int _sizeRowValY;

	FeatureMatrix *_FM;
	DContainer *_Y;
	DContainer *_X;
	DContainer *_valY;
	DContainer *_valX;
};

#endif
