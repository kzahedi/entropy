#include "FeatureMatrix.h"

FeatureMatrix::FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda){
		_valX= &eX;
		_valY= &eY;
		_X= &aX;
		_Y= &aY;
		_sizeColValY= _valY->columns();
		_sizeColValX= _valX->columns();
		_sizeRowValX= _valX->rows();
		_sizeRowValY= _valY->rows();
		_sizeX = _X->rows();
		_sizeY = _Y->rows();
		_FA=FeatureArray(valuelambda);
		_getMatrix(valuelambda);
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
		_sizeX = _X->rows();
		_sizeY = _Y->rows();
		_FA=FeatureArray(0);
		_getMatrix(0);
}
FeatureMatrix:: ~FeatureMatrix(){
		for(int i=0; i<_sizeColValX;i++){
			for(int j=0; j<_sizeColValY;j++){
				_FA[i][j].~Feature();
			}
		}
		delete _FA;

		//for(int i=0;i<_sizeRowValX;i++){
		//	for(int j=0;j<_sizeY;j++){
		//		_mat[i][j].clear();
		//	}
		//}
		//delete _mat;
}
// Indizes von Feature und von Lambda
double FeatureMatrix:: getFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY){
		assert(i<_sizeColValX && j<_sizeColValY);
		double lambda=_FA[i][j].getlambda(ilambdaX,ilambdaY);
		return lambda;
}
// Indizes von Feature und der Daten von X
double FeatureMatrix:: getFeatureArrayvalue(int i, int j,double ValX, double ValY){
		assert(i<_sizeColValX && j<_sizeColValY);
		double value=_FA[i][j].value(ValX,ValY);
		return value;
}
void FeatureMatrix::setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda){
		assert(i<_sizeColValX && j<_sizeColValY);
		_FA[i][j].setlambda(ilambdaX,ilambdaY,valuelambda);
}
int FeatureMatrix:: getFeatureArraydelta(int i, int j,int idelta, int jdelta, double ValX, double ValY){
		assert(i<_sizeColValX && j<_sizeColValY);
		assert(idelta < _sizeX && jdelta < _sizeY);
		int delta=_FA[i][j].delta(_X->get(idelta,0),_Y->get(jdelta,0), ValX,ValY);
		return delta;
}
vector<int> FeatureMatrix::getMatrixIndexX(int i, int j){
		assert(i<_sizeRowValX && j< _sizeRowValY );
		int k=0;
		vector<int> indX(0);
		while((*_Feati)(i,j,k)!=-2){
			indX.push_back((*_Feati)(i,j,k));
			k++;
		}
		return indX;
}
vector<int> FeatureMatrix:: getMatrixIndexY(int i, int j){
		assert(i<_sizeRowValX && j< _sizeY );
		vector<int> indY(0);
		int k=0;
		while((*_Featj)(i,j,k)!=-2){
			indY.push_back((*_Featj)(i,j,k));
			k++;
		}
		return indY;
}
vector<int> FeatureMatrix:: getMatrixIndexdX(int i,int j){
		assert(i<_sizeRowValX && j< _sizeY );
		vector<int> dindX(0);
		int k=0;
		while((*_Delti)(i,j,k)!=-2){
			dindX.push_back((*_Delti)(i,j,k));
			k++;
		}
		return dindX;
}
vector<int> FeatureMatrix:: getMatrixIndexdY(int i,int j){
		assert(i<_sizeRowValX && j< _sizeY );
		vector<int> dindY(0);
		int k=0;
		while((*_Deltj)(i,j,k)!=-2){
			dindY.push_back((*_Deltj)(i,j,k));
			k++;
		}
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
void FeatureMatrix:: _getMatrix(double valuelambda){
	_Feati= new SparseMatrix(-2);
	_Featj= new SparseMatrix(-2);
	_Delti= new SparseMatrix(-2);
	_Deltj= new SparseMatrix(-2);

	cout << "erstellen Matrix" << endl;
		for(int i=0;i<_sizeRowValX;i++){
			for(int j=0;j<_sizeY;j++){
				for(int varFeati=0;varFeati<_sizeColValX;varFeati++){
					for(int varFeatj=0;varFeatj<_sizeColValY;varFeatj++){
						if(_FA[varFeati][varFeatj].value((*_valX)(i,varFeati),(*_Y)(j,0))!=0){
							for(int deltai=0; deltai<_sizeX; deltai++ ){
								for(int deltaj=0; deltaj<_sizeY; deltaj++){
									if(_FA[varFeati][varFeatj].delta((*_X).get(deltai,0),(*_Y).get(deltaj,0),(*_valX).get(i, varFeati),(*_Y).get(j,0))!=-1){
										int k=0;
										while((*_Feati)(i,j,k)!=-2){k++;}
										(*_Feati)(i,j,k)=varFeati;
										(*_Featj)(i,j,k)=varFeatj;
										(*_Delti)(i,j,k)=deltai;
										(*_Deltj)(i,j,k)=deltaj;
									}
								}
							}
						}
					}
				}
			}
		}
		cout << " fertig " << endl;
}
