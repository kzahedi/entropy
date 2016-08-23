#include "FeatureMatrix.h"

FeatureMatrix::FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda){
		valX= &eX;
		valY= &eY;
		DContainer *X= &aX;
		DContainer *Y= &aY;
		_sizeColValY= (*valY).columns();
		_sizeColValX= (*valX).columns();
		_sizeRowValX= (*valX).columns();
		_sizeRowValY= (*valY).columns();
		_FA=FeatureArray(*valX,*valY,*X, *Y,valuelambda);
		getMatrix(*valX, *valY, valuelambda);
}
FeatureMatrix:: ~FeatureMatrix(){
		for(int i=0; i<_sizeColValX;i++){
			for(int j=0; j<_sizeColValY;j++){
				_FA[i][j].~Feature();
			}
		}
		delete _FA;

		for(int i=0;i<_sizeRowValX;i++){
			for(int j=0;j<_sizeRowValY;j++){
				_mat[i][j].clear();
			}
		}
		delete _mat;
}
// Indizes von Feature und von Lambda
double FeatureMatrix:: getFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY){
		assert(i<_sizeColValX && j<_sizeColValY);
		double lambda=_FA[i][j].getlambda(ilambdaX,ilambdaY);
		return lambda;
}
// Indizes von Feature und der Daten von X
double FeatureMatrix:: getFeatureArrayvalue(int i, int j,int RowValX, int RowValY){
		assert(i<_sizeColValX && j<_sizeColValY);
		double value=_FA[i][j].value((*valX)(RowValX,i),(*valY)(RowValY,j));
		return value;
}

void FeatureMatrix::setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double lambdavalue){
		assert(i<_sizeColValX && j<_sizeColValY);
		_FA[i][j].setlambda(ilambdaX,ilambdaY,lambdavalue);
}

vector<int> FeatureMatrix::getMatrixIndexX(int i, int j){
		assert(i<_sizeRowValX && j< _sizeRowValY );
		vector<int> indX = _mat[i][j][0];
		return indX;
}
vector<int> FeatureMatrix:: getMatrixIndexY(int i, int j){
		assert(i<_sizeRowValX && j< _sizeRowValY );
		vector<int> indY = _mat[i][j][1];
		return indY;
}


Feature** FeatureMatrix:: FeatureArray(DContainer &eX, DContainer &eY,DContainer &aX, DContainer &aY, double valuelambda){
		DContainer *valX= &eX;
		DContainer *valY= &eY;
		DContainer *X = &aX;
		DContainer *Y = &aY;
		Feature **FA;
		FA = new Feature*[_sizeColValX];
		for(int m=0; m< _sizeColValX; m++){
		  FA[m]= new Feature[_sizeColValY];
		}
		for(int m=0;m<_sizeColValX; m++){
		  for(int k=0;k<_sizeColValY;k++){
			  Feature *K=new Feature(*X,*Y,valuelambda);
			  FA[m][k]=*K;
		  }
		}
		return FA;
	}

void FeatureMatrix:: getMatrix(DContainer &eX, DContainer &eY,double valuelambda){
		DContainer *valX= &eX;
		DContainer *valY= &eY;
		vector<vector<int> > V(2,vector<int>(0));

		_mat = new vector<vector<int> >*[_sizeRowValX];
		for(int i=0;i<_sizeRowValX;i++){
			_mat[i]= new vector<vector<int> >[_sizeRowValY];
		}
		for(int i=0;i<_sizeRowValX;i++){
					for(int j=0;j<_sizeRowValY;j++){
								  _mat[i][j]= V;
					}
		}
		for(int i=0;i<_sizeRowValX;i++){
			for(int j=0;j<_sizeRowValY;j++){
				for(int varFeati=0;varFeati<_sizeColValX;varFeati++){
					for(int varFeatj=0;varFeatj<_sizeColValY;varFeatj++){
						if(_FA[varFeati][varFeatj].value((*valX)(i,varFeati),(*valY)(j,varFeatj))!=0){
							_mat[i][j][0].push_back(varFeati);
							_mat[i][j][1].push_back(varFeatj);
						}
					}
				}
			}
		}
}
