#include "SCGIS.h"

SCGIS::SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test){
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
    _FM= new InstanceMatrix(*_valX,*_valY,*_X,*_Y,lambdavalue);
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
    			_normaliser[i][j][k]=_sizeY;
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
    __scgis(maxit,konv,test);
}
double SCGIS:: getconv(int i){
	return _conv[i];
}
int SCGIS:: getsizeconv(){
	return _conv.size();
}
double SCGIS::getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY){
  assert(Feati<_sizeColValX && Featj<_sizeColValY);
  double lambda=_FM->getFeatureArraylambda(Feati,Featj, ilambdaX,ilambdaY);
  return lambda;
}
SCGIS:: ~SCGIS(){
	_FM->~InstanceMatrix();

	for(int i=0; i<_sizeColValX;i++){
		for(int j=0;j<_sizeColValY;j++){
			for(int k=0;k<_sizeX;k++){
				delete [] _observed[i][j][k];
			}
			delete [] _observed[i][j];
		}
		delete [] _observed[i];
	}
	delete [] _observed;

	for(int i=0;i< _sizeColValX;i++){
		for(int j=0;j<_sizeColValY;j++){
			for(int k=0;k<_sizeX;k++){
				delete [] _expected[i][j][k];
			}
			delete [] _expected[i][j];
		}
		delete [] _expected[i];
	}
	delete [] _expected;

	for(int i=0;i<_sizeColValX;i++){
		for(int j=0;j<_sizeColValY;j++){
			for(int k=0;k<_sizeRowValX;k++){
				delete [] _exponent[i][j][k];
			}
			delete [] _exponent[i][j];
		}
		delete [] _exponent[i];
	}
	delete [] _exponent;

	for(int i=0;i<_sizeColValX;i++){
		for(int j=0; j<_sizeColValY;j++){
			delete [] _normaliser[i][j];
		}
		delete [] _normaliser[i];
	}
	delete [] _normaliser;

	for(int i=0;i<_sizeX;i++){
		delete [] _delta[i];
	}
	delete [] _delta;

	_conv.clear();

}

double SCGIS::prop(int Feati,int Featj,double ValX,double ValY){
  double norm=0;
  double exponent=0;
  exponent= exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,ValY) );
  for(int yi=0;yi<_sizeY;yi++){
	    norm+= exp( (*_FM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
  }
  return exponent/norm;
}
double SCGIS::prop(int rowX,vector<vector<double> > Y, int rowY){
  double feat=0;
  double featnorm=0;
  double norm=0;
  double exponent=0;
  for(int Featx=0;Featx< _sizeColValX;Featx++){
	  for(int Featy=0;Featy< _sizeColValY;Featy++){
		  feat+=(*_FM).getFeatureArrayvalue(Featx,Featy,(*_valX)(rowX,Featx),Y[rowY][Featy]);
  	  }
    }
    exponent= exp(feat);
    for(int yi=0;yi<Y.size();yi++){
	  for(int Featx=0;Featx< _sizeColValX;Featx++){
		  for(int Featy=0;Featy< _sizeColValY;Featy++){
			  featnorm+= (*_FM).getFeatureArrayvalue(Featx,Featy,(*_valX)(rowX,Featx),Y[yi][Featy]);
		  }
	  }

	  norm+=exp(featnorm);
	  featnorm=0;
  }
  return exponent/norm;
}
double SCGIS:: prop(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY){
	  double feat=0;
	  double featnorm=0;
	  double norm=0;
	  double exponent=0;
	  for(int Featx=0;Featx< _sizeColValX;Featx++){
		  for(int Featy=0;Featy< _sizeColValY;Featy++){
			  feat+=(*_FM).getFeatureArrayvalue(Featx,Featy,X[rowX][Featx],Y[rowY][Featy]);
	  	  }
	    }
	    exponent= exp(feat);
	    for(int yi=0;yi<Y.size();yi++){
		  for(int Featx=0;Featx< _sizeColValX;Featx++){
			  for(int Featy=0;Featy< _sizeColValY;Featy++){
				  featnorm+= (*_FM).getFeatureArrayvalue(Featx,Featy,X[rowX][Featx],Y[yi][Featy]);
			  }
		  }

		  norm+=exp(featnorm);
		  featnorm=0;
	  }
	  return exponent/norm;
}

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
void SCGIS:: __scgis(int maxit, double konv, bool test){
	double l=1;
	int i=0;
while(i<maxit && fabs(l)>konv ){
	l=0;
	for(int Feati=0;Feati<_sizeColValX;Feati++){
		for(int Featj=0;Featj<_sizeColValY;Featj++){
			for(int delti=0;delti<_sizeX;delti++){
				for(int deltj=0;deltj<_sizeY; deltj++){
					_expected[Feati][Featj][delti][deltj]=0;
					for(int y=0;y<_sizeY;y++){
						for(int k=0; k<_FM->getInstanceMatrixX(Feati,Featj,delti,deltj).size();k++){
							if(_FM->getInstanceMatrixY(Feati,Featj,delti,deltj)[k]==y){
								int x=_FM->getInstanceMatrixX(Feati,Featj,delti,deltj)[k];
								if(fabs(_normaliser[Feati][Featj][x])>0.00000001 ){
									_expected[Feati][Featj][delti][deltj]+=exp(_exponent[Feati][Featj][x][y])/_normaliser[Feati][Featj][x];
								}
								else{cout << "norm " << _normaliser[Feati][Featj][x] << endl;}
							}
						}
					}
					double newl;
		            if(fabs(_expected[Feati][Featj][delti][deltj])<0.00000001){_expected[Feati][Featj][delti][deltj]=0.01;}
					if(fabs(_observed[Feati][Featj][delti][deltj])>0.00000001){
						_delta[delti][deltj]=log(_observed[Feati][Featj][delti][deltj]/_expected[Feati][Featj][delti][deltj]);
						newl= _FM->getFeatureArraylambda(Feati,Featj,delti,deltj)+_delta[delti][deltj];
					}
					else{newl=0; }
					l+=fabs( _observed[Feati][Featj][delti][deltj]-_expected[Feati][Featj][delti][deltj]);
					_FM->setFeatureArraylambda(Feati,Featj,delti,deltj,newl);
					for(int y=0;y<_sizeY;y++){
						for(int k=0; k<_FM->getInstanceMatrixX(Feati,Featj,delti,deltj).size();k++){
							if(_FM->getInstanceMatrixY(Feati,Featj,delti,deltj)[k]==y){
								int x=_FM->getInstanceMatrixX(Feati,Featj,delti,deltj)[k];
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
	 if(test){
		 _conv.push_back(l);
	 }
}
}
