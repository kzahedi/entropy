#ifndef __ITMATRIX_H__
#define __ITMATRIX_H__

#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"

using namespace std;

class ITMatrix
{
public:
	ITMatrix();
	ITMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue);
	virtual 	~ITMatrix();
 	double 		getFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY);
	double 		getFeatureArrayvalue(int i, int j,double ValX, double ValY);
	void 		setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda);
	int 		getFeatureArraydelta(int i, int j,int idelta, int jdelta, double ValX, double ValY);

protected:
	Feature** FeatureArray(double valuelambda);

	int _sizeColValX;
	int _sizeColValY;
	int _sizeRowValX;
	int _sizeRowValY;
	int _sizeX;
	int _sizeY;
	DContainer *_valX;
	DContainer *_valY;
	DContainer *_X;
	DContainer *_Y;
	Feature** _FA;
};

#endif
