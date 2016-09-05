#include "SCGIS.h"

SCGIS::SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,int maxit, double konv, double valuelambda){
    _valX= &eX;
    _valY= &eY;
    _X= &aX;
    _Y= &aY;
    _sizeX = (*_X).rows();
    _sizeY = (*_Y).rows();
    _sizeColValY= (*_valY).columns();
    _sizeColValX= (*_valX).columns();
    _sizeRowValX= (*_valX).rows();
    _sizeRowValY= (*_valY).rows();
    _FM=new InstanceMatrix(*_valX,*_valY,*_X,*_Y,valuelambda);
    _observed=__getobs();

    _exponent= new double***[_sizeColValX];
    for(int i=0;i<_sizeColValX; i++){
    	_exponent[i]=new double**[_sizeColValY];
    	for(int j=0;j<_sizeColValY;j++){
    		_exponent[i][j]=new double*[_sizeRowValX];
    		for(int xi=0;xi<_sizeRowValX;xi++){
    			_exponent[i][j][xi]=new double[_sizeY];
    			for(int y=0;y<_sizeY;y++){
    				_exponent[i][j][xi][y]=0;
    			}
    		}
    	}
    }

    _normaliser=new double**[_sizeColValX];
    for(int i=0; i< _sizeColValX;i++){
    	_normaliser[i]=new double*[_sizeColValY];
    	for(int j=0;j<_sizeColValY;j++){
    		_normaliser[i][j]=new double[_sizeRowValX];
    		for(int k=0;k<_sizeRowValX;k++){
    			_normaliser[i][j][k]=_sizeColValY;
    		}
    	}
    }
    _expected=new double***[_sizeColValX];
    for(int i=0;i<_sizeColValX;i++){
    	_expected[i]=new double**[_sizeColValY];
    	for(int j=0;j<_sizeColValY;j++){
    		_expected[i][j]=new double*[_sizeX];
    		for(int k=0;k<_sizeX;k++){
    			_expected[i][j][k]=new double[_sizeY];
    			for(int l=0;l<_sizeY;l++){
    				_expected[i][j][k][l]=0;
    			}
    		}
    	}
    }
    _delta= new double*[_sizeX];
    for(int i=0;i<_sizeX;i++){
    	_delta[i]=new double[_sizeY];
    	for(int j=0;j<_sizeY;j++){
    		_delta[i][j]=0;
    	}
    }
    __scgis(maxit,konv);
}
SCGIS:: ~SCGIS(){}

double**** SCGIS:: __getobs(){
    _observed = new double***[_sizeColValX];
    for(int i=0; i<_sizeColValX; i++){
      _observed[i]=new double**[_sizeColValY];
      for( int j=0;j< _sizeColValY;j++){
        _observed[i][j]=new double*[_sizeX];
        for(int k=0; k< _sizeX; k++){
          _observed[i][j][k]= new double[_sizeY];
          for(int l=0; l< _sizeY;l++){
            _observed[i][j][k][l]=0;
          }
        }
      }
    }

    //vector observed
    for(int i=0;i<_sizeRowValX;i++ ){
      for(int Feati=0; Feati< _sizeColValX;Feati++){
        for(int Featj=0;Featj< _sizeColValY; Featj++){
          for(int delti=0; delti< _sizeX; delti++){
            for(int deltj=0;deltj<_sizeY; deltj++){
              if(_FM->getFeatureArraydelta(Feati,Featj,delti,deltj,(*_valX)(i,Feati),(*_valY)(i,Featj))==1){
            		_observed[Feati][Featj][delti][deltj]++;
              }
            }
          }
        }
      }
    }
    return _observed;
}
void SCGIS:: __scgis(int maxit, double konv){
	int i=0;
while(i<maxit){
	for(int Feati=0;Feati<_sizeColValX;Feati++){
		for(int Featj=0;Featj<_sizeColValY;Featj++){
			for(int delti=0;delti<_sizeX;delti++){
				for(int deltj=0;deltj<_sizeY; deltj++){
					_expected[Feati][Featj][delti][deltj]=0;
					for(int y=0;y<_sizeY;y++){
						for(int k=0; k<_FM->getInstanceMatrixX(Feati,Featj,delti,deltj).size();k++){
							if(_FM->getInstanceMatrixY(Feati,Featj,delti,deltj)[k]==y){
								double x=_FM->getInstanceMatrixX(Feati,Featj,delti,deltj)[k];
								_expected[Feati][Featj][delti][deltj]+=exp(_exponent[Feati][Featj][x][y])/_normaliser[Feati][Featj][x];
							}
						}
					}
					_delta[delti][deltj]=log(_observed[Feati][Featj][delti][deltj]/_expected[Feati][Featj][delti][deltj]);
					double newl= _FM->getFeatureArraylambda(Feati,Featj,delti,deltj)+_delta[delti][deltj];
					_FM->setFeatureArraylambda(Feati,Featj,delti,deltj,newl);
					for(int y=0;y<_sizeY;y++){
						for(int k=0; k<_FM->getInstanceMatrixX(Feati,Featj,delti,deltj).size();k++){
							if(_FM->getInstanceMatrixY(Feati,Featj,delti,deltj)[k]==y){
								double x=_FM->getInstanceMatrixX(Feati,Featj,delti,deltj)[k];
								_normaliser[Feati][Featj][x]-=exp(_exponent[Feati][Featj][x][y]);
								_exponent[Feati][Featj][x][y]+=_delta[delti][deltj];
								_normaliser[Feati][Featj][x]+=exp(_exponent[Feati][Featj][x][y]);
							}
						}
					}
				}
			}
		}
	}
	i++;
}
}
