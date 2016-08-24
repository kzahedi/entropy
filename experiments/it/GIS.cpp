#include "GIS.h"

GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue) {
			DContainer *valX= &eX;
			DContainer *valY= &eY;
			DContainer *X= &aX;
			DContainer *Y= &aY;
			_sizeX = (*X).rows();
			_sizeY = (*Y).rows();
			_sizeColValY= (*valY).columns();
			_sizeColValX= (*valX).columns();
			_sizeRowValX= (*valX).rows();
			_sizeRowValY= (*valY).rows();
			FeatureMatrix *FM=new FeatureMatrix(*valX,*valY,*X,*Y,lambdavalue);
			gislambda(*FM);
}

GIS::~GIS() {

}

void GIS:: gislambda(FeatureMatrix &FM){
	double**** observed;
	observed = new double***[_sizeColValX];
	for(int i=0; i<_sizeColValX; i++){
		observed[i]=new double**[_sizeColValY];
		for( int j=0;j< _sizeColValY;j++){
			observed[i][j]=new double*[_sizeX];
			for(int k=0; k< _sizeX; k++){
				observed[i][j][k]= new double[_sizeY];
				for(int l=0; l< _sizeY;l++){
					observed[i][j][k][l]=0;
				}
			}
		}
	}

	//vector observed
	/*for(int i=0;i<_sizeColValX;i++ ){
		for(int j=0; j< _sizeColValY;j++){
			for(int k=0; k< FM.getMatrixIndexX(i,j).size();k++){
				observed[FM.getMatrixIndexX(i,j)[k]][FM.getMatrixIndexY(i,j)[k]]++;
				double m=observed[FM.getMatrixIndexX(i,j)[k]][FM.getMatrixIndexY(i,j)[k]];
			}
		}
	}
	//constant c
	int Featconst=0;
	int curr=0;
	for(int i=0; i< _sizeRowValX;i++){
		for(int j=0; j< _sizeRowValY;j++){
			for(int featxi=0; featxi < _sizeColValX; featxi++){
				for(int featyj=0; featyj < _sizeColValX; featyj++){
					for(int k=0; k< FM.getMatrixIndexX(i,j).size();k++){
						if(FM.getMatrixIndexX(i,j)[k]==featxi && FM.getMatrixIndexY(i,j)[k]==featyj){
							curr++;
						}
					}
				}
			}
		if(curr> Featconst) Featconst=curr;
		curr=0;
		}
	}
	*/







}

