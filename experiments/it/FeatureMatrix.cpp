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
		_FA=FeatureArray(valuelambda);
		getMatrix(valuelambda);
}
FeatureMatrix::FeatureMatrix(){
		_valX= new DContainer(0,0);
		_valY= new DContainer(0,0);
		_X= new DContainer(0,0);
		_Y= new DContainer(0,0);
		_sizeColValY= 0;
		_sizeColValX= 0;
		_sizeRowValX= 0;
		_sizeRowValY= 0;
		_sizeX = (*_X).rows();
		_sizeY = (*_Y).rows();
		_FA=FeatureArray(0);
		getMatrix(0);

}
FeatureMatrix:: ~FeatureMatrix(){
		for(int i=0; i<_sizeColValX;i++){
			for(int j=0; j<_sizeColValY;j++){
				_FA[i][j].~Feature();
			}
		}
		delete _FA;

		for(int i=0;i<_sizeRowValX;i++){
			for(int j=0;j<_sizeY;j++){
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
		double value=_FA[i][j].value((*_valX)(RowValX,i),(*_valY)(RowValY,j));
		return value;
}
double FeatureMatrix:: getFeatureArrayvalueval(int i, int j,double ValX, double ValY){
		assert(i<_sizeColValX && j<_sizeColValY);
		double value=_FA[i][j].value(ValX,ValY);
		return value;
}

void FeatureMatrix::setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda){
		assert(i<_sizeColValX && j<_sizeColValY);
		_FA[i][j].setlambda(ilambdaX,ilambdaY,valuelambda);
}
int FeatureMatrix:: getFeatureArraydelta(int i, int j,int idelta, int jdelta, int RowValX, int RowY){
		assert(i<_sizeColValX && j<_sizeColValY);
		assert(RowValX <_sizeRowValX && RowY <_sizeY);
		assert(idelta < _sizeX && jdelta < _sizeY);
		int delta=_FA[i][j].delta((*_X).get(idelta,0),(*_Y).get(jdelta,0), (*_valX)(RowValX,i),(*_Y)(RowY,0));
		return delta;
}
int FeatureMatrix:: getFeatureArraydeltaval(int i, int j,int idelta, int jdelta, double ValX, double ValY){
		assert(i<_sizeColValX && j<_sizeColValY);
		assert(idelta < _sizeX && jdelta < _sizeY);
		int delta=_FA[i][j].delta((*_X).get(idelta,0),(*_Y).get(jdelta,0), ValX,ValY);
		return delta;
}
vector<int> FeatureMatrix::getMatrixIndexX(int i, int j){
		assert(i<_sizeRowValX && j< _sizeRowValY );
		vector<int> indX = _mat[i][j][0];
		return indX;
}
vector<int> FeatureMatrix:: getMatrixIndexY(int i, int j){
		assert(i<_sizeRowValX && j< _sizeY );
		vector<int> indY = _mat[i][j][1];
		return indY;
}
vector<int> FeatureMatrix:: getMatrixIndexdX(int i,int j){
		assert(i<_sizeRowValX && j< _sizeY );
		vector<int> dindX = _mat[i][j][2];
		return dindX;
}
vector<int> FeatureMatrix:: getMatrixIndexdY(int i,int j){
		assert(i<_sizeRowValX && j< _sizeY );
		vector<int> dindY = _mat[i][j][3];
		return dindY;
}

Feature** FeatureMatrix:: FeatureArray( double valuelambda){
		Feature **FA;
		FA = new Feature*[_sizeColValX];
		for(int m=0; m< _sizeColValX; m++){
		  FA[m]= new Feature[_sizeColValY];
		  for(int k=0;k<_sizeColValY;k++){
			  Feature *K=new Feature(*_X,*_Y,valuelambda);
			  FA[m][k]=*K;
		  }
		}
		return FA;
}

void FeatureMatrix:: getMatrix(double valuelambda){
		vector<vector<int> > V(4,vector<int>(0));

		_mat = new vector<vector<int> >*[_sizeRowValX];
		for(int i=0;i<_sizeRowValX;i++){
			_mat[i]= new vector<vector<int> >[_sizeY];
			for(int j=0;j<_sizeY;j++){
				_mat[i][j]= V;
			}
		}

		for(int i=0;i<_sizeRowValX;i++){
			for(int j=0;j<_sizeY;j++){
				for(int varFeati=0;varFeati<_sizeColValX;varFeati++){
					for(int varFeatj=0;varFeatj<_sizeColValY;varFeatj++){
						if(_FA[varFeati][varFeatj].value((*_valX)(i,varFeati),(*_Y)(j,0))!=0){
							for(int deltai=0; deltai<_sizeX; deltai++ ){
								for(int deltaj=0; deltaj<_sizeY; deltaj++){
									if(_FA[varFeati][varFeatj].delta((*_X).get(deltai,0),(*_Y).get(deltaj,0),(*_valX).get(i, varFeati),(*_valY).get(j,0))!=0 ){
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
