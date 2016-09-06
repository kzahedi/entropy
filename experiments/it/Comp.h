#ifndef __COMP_H__
#define __COMP_H__

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

class Comp{

public:
	Comp(GIS &exact, DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,int maxit, double konv);
	~Comp();
	vector<double> comptimeKL(int maxit, double konv, vector<double> y);
	vector<double> KL(vector<double> y );
	//uebersicht

private:
	GIS 	*_exact;
	GIS 	*_gisTest;
	SCGIS 	*_scgisTest;
	DContainer 		*_Y;
	DContainer 		*_X;
	DContainer 		*_valY;
	DContainer 		*_valX;
	vector<double>	_timediff;

	vector<vector<double> > getY();
	double 	__prop(GIS &test, int RowX, vector<double> Y);
	double 	__prop(SCGIS &test, int RowX,vector<double> Y);
	void 	__comptime(int maxit, double konv);
};

#endif
