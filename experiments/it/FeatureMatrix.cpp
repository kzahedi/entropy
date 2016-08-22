#include "FeatureMatrix.h"
#include "Feature.h"


FeatureMatrix::FeatureMatrix(DContainer &aX, DContainer &aY, DContainer &eX, DContainer &eY, double valuelambda){
		DContainer *X= &aX;
		DContainer *Y= &aY;
		FA=getFeatures(*X, *Y,valuelambda);
	}

Feature** FeatureMatrix:: getFeatures(DContainer &aX, DContainer &aY, double valuelambda){
		DContainer *X= &aX;
		DContainer *Y= &aY;
		_sizeY= (*Y).rows();
		_sizeX= (*X).rows();
		Feature **FA;
		FA = new Feature*[_sizeX];
		for(int m=0; m< _sizeX; m++){
		  FA[m]= new Feature[_sizeY];
		}
		for(int m=0;m<_sizeX; m++){
		  for(int k=0;k<_sizeY;k++){
			  Feature *K=new Feature(*X,*Y,valuelambda);
			  FA[m][k]=*K;
		  }
		}
		return FA;
	}
 //#include <tuple>
// typedef std::tuple<int,int> Tuple;
// Tuple** matrix;
