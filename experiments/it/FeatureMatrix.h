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
	Feature** FA;
	FeatureMatrix(DContainer &aX, DContainer &aY,double la);

private:
	Feature** getArray(DContainer &aX, DContainer &aY, double la);
	int _sizeX;
	int _sizeY;

};

#endif
