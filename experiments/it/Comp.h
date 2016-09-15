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
#include "GISgp.h"
#include "SCGISgp.h"

using namespace std;

class Comp{

public:
	Comp(int ColX,int RowX,int ColValY,  vector<double> lambda,DContainer &aX, DContainer &aY,int maxit, double konv,bool time, int seconds);
	~Comp();
	vector<double> KL(vector< vector<double> > y);
	void comparison();

private:
	GIS 	*_exact;
	GIS 	*_gisTest;
	GISgp	*_gisgpTest;
	SCGIS 	*_scgisTest;
	SCGISgp	*_scgisgpTest;
	DContainer 		*_Y;
	DContainer 		*_X;
	DContainer 		*_valY;
	DContainer 		*_valX;
	int 			_sizeColValY;
	int				_sizeColValX;
	int 			_index;
	vector<double>	_timediff;
	vector<vector<double> > _alphY;
	vector<vector<double> > _alphX;
	bool			_timetest;

	vector<vector<double> > __getalph(bool valX);

	void			__getValY(int ColY,int RowX);
	void		 	__getValX(int ColX,int RowX);
	void 			__setlambda(IContainer &indizes, DContainer &values);
	void 			__setlambdarand(vector<double> lambdaval);
	void 			__comptimegiss(int maxit);
	void 			__comptime(int maxit, double konv);
	void 			__comptime(int maxit, double konv,int seconds);
	void    		__fill(vector<double> fill, int i, vector<vector<double> > &Z,bool valX,int rowsAlph,int colval);
};

#endif
