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
	GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double la);
	double gis(int x,int y);
	~GIS();
private:
	void __gis();
	double**** __getobs();
	int _sizeX;
	int _sizeY;
	int _sizeColValX;
	int _sizeColValY;
	int _sizeRowValX;
	int _sizeRowValY;
	FeatureMatrix *_FM;
	DContainer *_valY;
	DContainer *_valX;
	double** __getFeatconst();
	void __getexp(double**** &expect, double*** &exponent,double** &normaliser);
};

#endif
