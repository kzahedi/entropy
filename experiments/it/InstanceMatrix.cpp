#include "InstanceMatrix.h"

InstanceMatrix:: InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda){
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
	_getMatrix(valuelambda);
}
InstanceMatrix:: ~InstanceMatrix(){
	// matrix
	for(int i=0; i<_sizeColValX;i++){
		for(int j=0; j<_sizeColValY;j++){
			_FA[i][j].~Feature();
		}
	}
	delete _FA;
}
Feature** InstanceMatrix :: FeatureArray( double valuelambda){
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
// Indizes von Feature und von Lambda
double InstanceMatrix:: getFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY){
		assert(i<_sizeColValX && j<_sizeColValY);
		double lambda=_FA[i][j].getlambda(ilambdaX,ilambdaY);
		return lambda;
}
// Indizes von Feature und der Daten von X
double InstanceMatrix:: getFeatureArrayvalue(int i, int j,double ValX, double ValY){
		assert(i<_sizeColValX && j<_sizeColValY);
		double value=_FA[i][j].value(ValX,ValY);
		return value;
}

void InstanceMatrix::setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda){
		assert(i<_sizeColValX && j<_sizeColValY);
		_FA[i][j].setlambda(ilambdaX,ilambdaY,valuelambda);
}
int InstanceMatrix:: getFeatureArraydelta(int i, int j,int idelta, int jdelta, double ValX, double ValY){
		assert(i<_sizeColValX && j<_sizeColValY);
		assert(idelta < _sizeX && jdelta < _sizeY);
		int delta=_FA[i][j].delta((*_X).get(idelta,0),(*_Y).get(jdelta,0), ValX,ValY);
		return delta;
}
void InstanceMatrix::_getMatrix(double valuelambda){
		vector<vector<int> > V(2,vector<int>(0));
		_mat = new vector<vector<int> >***[_sizeColValX];
		for(int i=0;i<_sizeColValX;i++){
			_mat[i]= new vector<vector<int> >**[_sizeColValY];
			for(int j=0;j<_sizeColValY;j++){
				_mat[i][j]= new vector<vector<int> >*[_sizeX];
				for(int k=0;k<_sizeX;k++){
					_mat[i][j][k]= new vector<vector<int> >[_sizeY];
					for(int l=0;l<_sizeY;l++){
						_mat[i][j][k][l]=V;
					}
				}
			}
		}
		for(int Feati=0;Feati<_sizeColValX;Feati++){
			for(int Featj=0;Featj<_sizeColValY;Featj++){
				for(int delti=0; delti<_sizeX;delti++){
					for(int deltj=0; deltj < _sizeY; deltj++){
						for(int xi=0; xi< _sizeRowValX; xi++){
							for(int y=0; y<_sizeY; y++){
								if(_FA[Feati][Featj].delta(_X->get(delti,0),_Y->get(deltj,0),_valX->get(xi,Feati),_Y->get(y,0))==1){
									_mat[Feati][Featj][delti][deltj][0].push_back(xi);
									_mat[Feati][Featj][delti][deltj][1].push_back(y);
								}
							}
						}
					}
				}
			}
		}

}
