#include "InstanceMatrix.h"

InstanceMatrix:: InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda)
:ITMatrix(eX,eY,aX,aY,valuelambda){
	_getMatrix(valuelambda);
}
InstanceMatrix:: ~InstanceMatrix(){
	for(int i=0; i<_sizeColValX;i++){
		for(int j=0; j<_sizeColValY;j++){
			_FA[i][j].~Feature();
		}
	}
	delete _FA;
	for(int i=0;i<_sizeColValX;i++){
		for(int j=0;j<_sizeColValY;j++){
			for(int k=0;k<_sizeX;k++){
				for(int l=0;l<_sizeY;l++){
					_mat[i][j][k][l].clear();
				}
			}
		}
	}
	delete _mat;
}
vector<int> InstanceMatrix::getInstanceMatrixX(int Feati,int Featj,int delti,int deltj){
	assert(Feati<_sizeColValX && Featj<_sizeColValY);
	assert(delti<_sizeX && deltj< _sizeY);
	return _mat[Feati][Featj][delti][deltj][0];

}
vector<int> InstanceMatrix::getInstanceMatrixY(int Feati,int Featj,int delti,int deltj){
	assert(Feati<_sizeColValX && Featj<_sizeColValY);
	assert(delti<_sizeX && deltj< _sizeY);
	return _mat[Feati][Featj][delti][deltj][1];
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
