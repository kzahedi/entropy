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

	FeatureMatrix(DContainer &aX, DContainer &aY, DContainer &eX, DContainer &eY,double la);
	~FeatureMatrix();
	double getFeatureArraylambda(int i, int j,int k, int l);
	double getFeatureArrayvalue(int i, int j, int k, int l);
	void setFeatureArraylambda(int i, int j,int k, int l,double lambdavalue);
	vector<int> getMatrixIndexX(int i, int j);
	vector<int> getMatrixIndexY(int i, int j);

private:
	Feature** FeatureArray(DContainer &eX, DContainer &eY,DContainer &aX, DContainer &aY, double la);
	void getMatrix(DContainer &eX, DContainer &eY,double la);

	int _sizeColValX;
	int _sizeColValY;
	int _sizeRowValX;
	int _sizeRowValY;
	DContainer *valX;
	DContainer *valY;
	Feature** _FA;
	vector<vector<int> > **_mat;
};

#endif
