#ifndef __FEATUREMATRIX_H__
#define __FEATUREMATRIX_H__

#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

class FeatureMatrix
{
public:
		Feature** FA;

	FeatureMatrix(DContainer &aX, DContainer &aY,double la){
		DContainer *X= &aX;
		DContainer *Y= &aY;
		FA=getArray(*X, *Y,la);
	}

private:
	Feature** getArray(DContainer &aX, DContainer &aY, double la){
		DContainer *X= &aX;
		DContainer *Y= &aY;
		sizeY= (*Y).rows();
		sizeX= (*X).rows();
		Feature **FA;
		FA = new Feature*[sizeX];
		for(int m=0; m< sizeX; m++){
		  FA[m]= new Feature[sizeY];
		}
		for(int m=0;m<sizeX; m++){
		  for(int k=0;k<sizeY;k++){
			  Feature *K=new Feature(*X,*Y,la);
			  FA[m][k]=*K;
		  }
		}
		return FA;
	}
	int sizeX;
	int sizeY;

};

#endif
