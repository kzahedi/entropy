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
	for(int i=0;i<_sizeColValX;i++ ){
		for(int j=0; j< _sizeColValY;j++){
			for(int k=0; k< FM.getMatrixIndexX(i,j).size();k++){
				observed[FM.getMatrixIndexX(i,j)[k]][FM.getMatrixIndexY(i,j)[k]][FM.getMatrixIndexdX(i,j)[k]][FM.getMatrixIndexdY(i,j)[k]]++;
			}
		}
	}
	//constant c for delta
	int** Featconst;
	Featconst = new int*[_sizeColValX];
	for(int i=0; i< _sizeColValX;i++){
		Featconst[i]=new int[_sizeColValY];
		for(int j=0; j< _sizeColValY; j++){
			Featconst[i][j]=0;
		}
	}
	int curr=0;
	for(int delti=0; delti< _sizeColValX; delti++){
		for(int deltj=0; deltj< _sizeColValY; deltj++){

			for(int i=0; i< _sizeRowValX;i++){
				for(int j=0; j< _sizeRowValY;j++){

					for(int deltxi=0; deltxi < _sizeX; deltxi++){
						for(int deltyj=0; deltyj < _sizeY; deltyj++){
							for(int k=0; k< FM.getMatrixIndexX(i,j).size();k++){
								if(FM.getMatrixIndexdX(i,j)[k]==deltxi && FM.getMatrixIndexdY(i,j)[k]==deltyj){
									curr++;
								}
							}
						}
					}
					if(curr> Featconst[delti][deltj]) Featconst[delti][deltj]=curr;
					curr=0;
				}
			}
		}
	}








}

