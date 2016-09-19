#ifndef __INSTANCEMATRIX_H__
#define __INSTANCEMATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>
#include <entropy++/Matrix.h>
#include "Feature.h"

using namespace std;

class InstanceMatrix{

public:
	//InstanceMatrix();
	InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda);
	~InstanceMatrix();
	vector<int> getInstanceMatrixX(int Feati,int Featj,int delti,int deltj);
	vector<int> getInstanceMatrixY(int Feati,int Featj,int delti,int deltj);
	double 	getFeatureArraylambda(int Feati,int Featj,int ilambdaX,int ilambdaY);
	double 	getFeatureArrayvalue(int Feati,int Featj,double ValX,double ValY);
	void 	setFeatureArraylambda(int Feati,int Featj,int ilambdaX,int ilambdaY,double valuelambda);
	int 	getFeatureArraydelta(int Feati,int Featj,int idelta,int jdelta, double ValX, double ValY);

private:
	Feature** FeatureArray(double valuelambda);
	void _getMatrix(double valuelambda);

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
	vector<vector<int> > ****_mat;

};
#endif
