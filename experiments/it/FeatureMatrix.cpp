#include "FeatureMatrix.h"

FeatureMatrix::FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda){
		_valX= &eX;
		_valY= &eY;
		_X= &aX;
		_Y= &aY;
		_sizeColValY= (*_valY).columns();
		_sizeColValX= (*_valX).columns();
		_sizeRowValX= (*_valX).rows();
		_sizeRowValY= (*_valY).rows();
		_sizeX = (*_X).rows();
		_sizeY = (*_Y).rows();
		_FA=FeatureArray(*_valX,*_valY,*_X, *_Y,valuelambda);
		getMatrix(*_valX, *_valY, valuelambda);
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
		_X=NULL;
		_Y=NULL;
		_valX=NULL;
		_valY=NULL;
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
		double value=_FA[i][j].value((*_valX)(RowValX,i),(*_valY)(RowValY,j));
		return value;
}
double FeatureMatrix:: getFeatureArrayvalueforval(int i, int j,int x, int y){
		assert(i<_sizeColValX && j<_sizeColValY);
		double value=_FA[i][j].value(x,y);
		return value;
}

void FeatureMatrix::setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda){
		assert(i<_sizeColValX && j<_sizeColValY);
		_FA[i][j].setlambda(ilambdaX,ilambdaY,valuelambda);
}
int FeatureMatrix:: getFeatureArraydelta(int i, int j,int idelta, int jdelta, int RowValX, int RowValY){
		assert(i<_sizeColValX && j<_sizeColValY);
		assert(RowValX <_sizeRowValX && RowValY <_sizeRowValY);
		assert(idelta < _sizeX && jdelta << _sizeY);
		int delta=_FA[i][j].delta((*_X).get(idelta,0),(*_Y).get(jdelta,0), (*_valX)(RowValX,i),(*_valY)(RowValY,j));
		return delta;
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
vector<int> FeatureMatrix:: getMatrixIndexdX(int i,int j){
		assert(i<_sizeRowValX && j< _sizeRowValY );
		vector<int> dindX = _mat[i][j][2];
		return dindX;
}
vector<int> FeatureMatrix:: getMatrixIndexdY(int i,int j){
		assert(i<_sizeRowValX && j< _sizeRowValY );
		vector<int> dindY = _mat[i][j][3];
		return dindY;
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
		vector<vector<int> > V(4,vector<int>(0));

		_mat = new vector<vector<int> >*[_sizeRowValX];
		for(int i=0;i<_sizeRowValX;i++){
			_mat[i]= new vector<vector<int> >[_sizeRowValY];
			for(int j=0;j<_sizeRowValY;j++){
				_mat[i][j]= V;
			}
		}

		for(int i=0;i<_sizeRowValX;i++){
			for(int j=0;j<_sizeRowValY;j++){
				for(int varFeati=0;varFeati<_sizeColValX;varFeati++){
					for(int varFeatj=0;varFeatj<_sizeColValY;varFeatj++){
						if(_FA[varFeati][varFeatj].value((*valX)(i,varFeati),(*valY)(j,varFeatj))!=0){
							for(int deltai=0; deltai<_sizeX; deltai++ ){
								for(int deltaj=0; deltaj<_sizeY; deltaj++){
									if(_FA[varFeati][varFeatj].delta((*_X).get(deltai,0),(*_Y).get(deltaj,0),(*_valX).get(i, varFeati),(*_valY).get(j,varFeatj))!=0 ){
										_mat[i][j][0].push_back(varFeati);
										_mat[i][j][1].push_back(varFeatj);
										_mat[i][j][2].push_back(deltai);
										_mat[i][j][3].push_back(deltaj);
									}
								}
							}
						}
					}
				}
			}
		}
}
