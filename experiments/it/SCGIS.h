#ifndef __SCGIS_H__
#define __SCGIS_H__

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

class SCGIS{

public:
	SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,int maxit, double konv, double valuelambda);
	~SCGIS();

private:
	void			__scgis(int maxit, double konv);
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
	double**		_delta;

	InstanceMatrix 	*_FM;
	DContainer 		*_Y;
	DContainer 		*_X;
	DContainer 		*_valY;
	DContainer 		*_valX;

};
#endif
