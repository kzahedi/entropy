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

using namespace std;

class InstanceMatrix{

public:
	//InstanceMatrix();
	InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda);
	~InstanceMatrix();
	// Matrix indizes Abfragen vier eingaben, vektorx als ausgabe
	// Matrix abfragen, vier eingaben, vektor y als ausgabe
	double 	getFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY);
	double 	getFeatureArrayvalue(int i, int j,double ValX, double ValY);
	void 	setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda);
	int 	getFeatureArraydelta(int i, int j,int idelta, int jdelta, double ValX, double ValY);

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
