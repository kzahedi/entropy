#ifndef __FEATUREMATRIX_H__
#define __FEATUREMATRIX_H__

#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"

using namespace std;

class FeatureMatrix
{
public:
	Feature** FA;   																				//ueber get
	FeatureMatrix(DContainer &aX, DContainer &aY, DContainer &eX, DContainer &eY,double la);
private:
	Feature** getFeatures(DContainer &aX, DContainer &aY, double la);
	void getMatrix(DContainer &eX, DContainer &eY,double la);
	int _sizeX;
	int _sizeY;


};

#endif
