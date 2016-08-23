#include "FeatureMatrix.h"

FeatureMatrix::FeatureMatrix(DContainer &aX, DContainer &aY, DContainer &eX, DContainer &eY, double valuelambda){
		DContainer *valX= &eX;
		DContainer *valY= &eY;
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

void FeatureMatrix:: getMatrix(DContainer &eX, DContainer &eY,double la){
		DContainer *valX= &eX;
		DContainer *valY= &eY;
		int sizeValX=(*valX).rows();
		int sizeValY=(*valY).rows();
		int sizecolX=(*valX).columns();
		int sizecolY=(*valY).columns();
		vector<vector<int> > V(2,vector<int>(0));
		vector<vector<int> >mat[sizeValX][sizeValY];
				for(int i=0;i<sizeValX;i++){
					  for(int j=0;j<sizeValY;j++){
						  mat[i][j]= V;
					  }
				}
		for(int i=0;i<sizeValX;i++){
			for(int j=0;j<sizeValY;j++){
				for(int varFeati=0;varFeati<sizecolX;varFeati++){
					for(int varFeatj=0;varFeatj<sizecolY;varFeatj++){
						if(FA[varFeati][varFeatj].value((*valX)(i,varFeati),(*valY)(j,varFeatj))!=0){
							mat[i][j][0].push_back(varFeati);
							mat[i][j][1].push_back(varFeatj);
						}
					}
				}
			}
		}
}
