#ifndef __IT_H__
#define __IT_H__

#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_MI.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>
#include <entropy++/sparse/state/MC_MI.h>
#include <entropy++/Matrix.h>
#include <entropy++/SparseMatrix.h>

#include <time.h>
#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"
#include "ITMatrix.h"
#include "FeatureMatrix.h"
#include "InstanceMatrix.h"

using namespace std;

class IT{

public:
	IT(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,bool GIS);
	IT(int ColValY, DContainer &eX,DContainer &aX, DContainer &aY);
	double 	prop(int Feati,int Featj,double ValX,double ValY);
	double 	prop(int rowX,vector<vector<double> > Y, int rowY);
	double 	prop(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY);
	double 	propm(vector<vector<double> > X,int rowX,vector<vector<double> > Y);
	double 	propm(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY);
	double	getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY);

protected:
	double**** 		__getobs();
	int 			_sizeX;
	int 			_sizeY;
	int 			_sizeColValX;
	int 			_sizeColValY;
	int 			_sizeRowValX;
	int 			_sizeRowValY;
	double**** 		_observed;
	bool			_gis;

	DContainer 		*_Y;
	DContainer 		*_X;
	DContainer 		*_valY;
	DContainer 		*_valX;
	FeatureMatrix	*_FM;
	InstanceMatrix  *_IM;
};
#endif
