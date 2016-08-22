#include "Feature.h"

Feature::Feature() {
	X = new DContainer(0, 0);
	Y = new DContainer(0, 0);
	sizeY = (*Y).rows();
	_sizeX = (*X).rows();
	lambda = new double*[_sizeX];
	for (int m = 0; m < _sizeX; m++)
		lambda[m] = new double[sizeY];
}

double Feature::getlambda(int i, int j) {
	assert(i < _sizeX && j < sizeY);
	return lambda[i][j];
}

void Feature::setlambda(int i, int j, double newvalue) {
	assert(i < _sizeX && j < sizeY);
	lambda[i][j] = newvalue;
}


#include <tuple>

typedef std::tuple<int,char> Tuple;
Tuple** matrix;

