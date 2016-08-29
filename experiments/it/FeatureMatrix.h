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
	FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue);
	FeatureMatrix();
	~FeatureMatrix();
	double getFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY);
	double getFeatureArrayvalue(int i, int j,int RowValX, int RowValY);
	void setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda);
	int getFeatureArraydelta(int i, int j,int idelta, int jdelta, int RowValX, int RowValY);
	vector<int> getMatrixIndexX(int i, int j);
	vector<int> getMatrixIndexY(int i, int j);
	vector<int> getMatrixIndexdX(int i,int j);
	vector<int> getMatrixIndexdY(int i,int j);

private:
	Feature** FeatureArray(double la);
	void getMatrix(double la);

	int _sizeColValX;
	int _sizeColValY;
	int _sizeRowValX;
	int _sizeRowValY;
	int _sizeX;
	int _sizeY;
	DContainer *_valX;
	DContainer *_valY;
	DContainer *_X;
	DContainer *_Y;
	Feature** _FA;
	vector<vector<int> > **_mat;
};

#endif
