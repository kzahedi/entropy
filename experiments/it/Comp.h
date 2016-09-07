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
	vector<double> KL(vector< vector<double> > y, int RowY );
	void comparison(int RowY);

private:
	GIS 	*_exact;
	GIS 	*_gisTest;
	SCGIS 	*_scgisTest;
	DContainer 		*_Y;
	DContainer 		*_X;
	DContainer 		*_valY;
	DContainer 		*_valX;
	int 			_sizeColValY;
	int				_sizeColValX;
	vector<double>	_timediff;
	vector<vector<double> > _alphY;
	vector<vector<double> > _alphX;

	vector<vector<double> > __getY();
	vector<vector<double> > __getX();

	void 	__comptime(int maxit, double konv);
	void    __fill(vector<double> fill, int i, vector<vector<double> > &Y);
	void    __fillx(vector<double> fill, int i, vector<vector<double> > &Y);
};

#endif
