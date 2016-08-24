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
	~GIS();
private:
	void gislambda(FeatureMatrix &FM);
	int _sizeColValX;
	int _sizeColValY;
	int _sizeRowValX;
	int _sizeRowValY;
};

#endif
