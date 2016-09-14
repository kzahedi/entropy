#include "GIS.h"


GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test) {
      _valX= &eX;
      _valY= &eY;
      _X= &aX;
      _Y= &aY;
      _sizeX = _X->rows();
      _sizeY = _Y->rows();
      _sizeColValY= _valY->columns();
      _sizeColValX= _valX->columns();
      _sizeRowValX= _valX->rows();
      _sizeRowValY= _valY->rows();
      assert(_valX->rows()==_valY->rows());
      _FM=new FeatureMatrix(*_valX,*_valY,*_X,*_Y,lambdavalue);


	  _exponent = new double**[_sizeColValX];
	  for(int i=0; i<_sizeColValX; i++){
	    _exponent[i]=new double*[_sizeColValY];
	    for(int j=0;j<_sizeColValY;j++){
	      _exponent[i][j]=new double[_sizeRowValY];
	      }
	    }

	  _normaliser= new double*[_sizeColValX];
	  for(int i=0; i< _sizeColValX; i++){
	    _normaliser[i]=new double[_sizeColValY];
	  }

	  _expected = new double***[_sizeColValX];
	  for(int i=0; i<_sizeColValX; i++){
	    _expected[i]=new double**[_sizeColValY];
	    for( int j=0;j< _sizeColValY;j++){
	      _expected[i][j]=new double*[_sizeX];
	      for(int k=0; k< _sizeX; k++){
	        _expected[i][j][k]= new double[_sizeY];
	        for(int l=0; l< _sizeY;l++){
	          _expected[i][j][k][l]=0;
	        }
	      }
	    }
	  }
	  _observed=__getobs();
	  __gis(maxit, konv,true);
}

// ColValY - Anzahl der Y-Knoten
GIS::GIS(int ColValY, DContainer &eX,DContainer &aX, DContainer &aY){
    _X= &aX;
    _Y= &aY;
    _sizeX=_X->rows();
    _sizeY=_Y->rows();
    _valX= &eX;
    _sizeColValY= ColValY;
    _sizeColValX= (*_valX).columns();
    _sizeRowValX= (*_valX).rows();
    _valY= new DContainer(_sizeRowValX,ColValY);
    _sizeRowValY= 0;
    _FM=new FeatureMatrix(*_valX,*_valY,*_X,*_Y,0);

    _exponent = NULL;
    _normaliser= NULL;
    _expected = NULL;
    _observed = NULL;

}

GIS::~GIS(){
	if(_observed!=NULL){
	for(int i=0;i<_sizeColValX;i++){
	   for(int j=0;j<_sizeColValY;j++){
	      for(int k=0;k<_sizeX;k++){
	         delete [] _observed[i][j][k];
		  }
		  delete [] _observed[i][j];
	   }
	   delete [] _observed[i];
	}
	delete [] _observed;
	}
	if(_expected != NULL){
	   for(int i=0;i<_sizeColValX;i++){
	      for(int j=0;j<_sizeColValY;j++){
		     for(int k=0;k<_sizeX;k++){
			    delete [] _expected[i][j][k];
			 }
			 delete [] _expected[i][j];
		  }
		  delete[] _expected[i];
	   }
	   delete[] _expected;
	}
	if(_normaliser != NULL){
		for(int m=0;m<_sizeColValX;m++){
			delete [] _normaliser[m];
		}
		delete [] _normaliser;
	}
	if(_exponent != NULL){
		for(int i=0;i<_sizeColValX;i++){
			for(int j=0;j<_sizeColValY;j++){
				delete [] _exponent[i][j];
			}
			delete [] _exponent[i];
		}
		delete [] _exponent;
	}
	_FM->~FeatureMatrix();
	_conv.clear();
}
double GIS::prop(int Feati,int Featj,double ValX,double ValY){
  double norm=0;
  double exponent=0;
  exponent= exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,ValY) );
  for(int yi=0;yi<_sizeY;yi++){
	    norm+= exp( (*_FM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
  }
  return exponent/norm;
}
double GIS::prop(int rowX,vector<vector<double> > Y, int rowY){
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
double GIS:: prop(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY){
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
double GIS::	propm(vector<vector<double> > X,int rowX,vector<vector<double> > Y){
	double z=0;
	double feat=0;
	for(int y=0;y<Y.size();y++){
		  for(int Featx=0;Featx< _sizeColValX;Featx++){
			  for(int Featy=0;Featy< _sizeColValY;Featy++){
				  feat+=(*_FM).getFeatureArrayvalue(Featx,Featy,X[rowX][Featx],Y[y][Featy]);
		  	  }
		    }
		  z+=exp(feat);
		  feat=0;
	}
	double n;
	double nexp=0;
	for(int x=0;x<X.size();x++){
		for(int y=0;y<Y.size();y++){
			  for(int Featx=0;Featx< _sizeColValX;Featx++){
				  for(int Featy=0;Featy< _sizeColValY;Featy++){
					  feat+=(*_FM).getFeatureArrayvalue(Featx,Featy,X[x][Featx],Y[y][Featy]);
			  	  }
			  }
			  nexp +=exp(feat);
			  feat=0;
		}
		n+=nexp;
		nexp=0;
	}
	return z/n;
}
double GIS:: propm(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY){
	  double feat=0;
	  double featnorm=0;
	  double norm=0;
	  double n=0;
	  double exponent=0;
	  for(int Featx=0;Featx< _sizeColValX;Featx++){
		  for(int Featy=0;Featy< _sizeColValY;Featy++){
			  feat+=(*_FM).getFeatureArrayvalue(Featx,Featy,X[rowX][Featx],Y[rowY][Featy]);
	  	  }
	    }
	   exponent= exp(feat);
	   for(int x=0;x<X.size();x++){
	    for(int yi=0;yi<Y.size();yi++){
		  for(int Featx=0;Featx< _sizeColValX;Featx++){
			  for(int Featy=0;Featy< _sizeColValY;Featy++){
				  featnorm+= (*_FM).getFeatureArrayvalue(Featx,Featy,X[rowX][Featx],Y[yi][Featy]);
			  }
		  }

		  norm+=exp(featnorm);
		  featnorm=0;
	    }
	    n+=norm;
	    norm=0;
	   }
	  return exponent/n;
}
double GIS:: getconv(int i){
	return _conv[i];
}
int GIS:: getsizeconv(){
	return _conv.size();
}
void GIS::setFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY,double valuelambda){
  assert(Feati<_sizeColValX && Featj<_sizeColValY);
  _FM->setFeatureArraylambda(Feati,Featj,ilambdaX,ilambdaY,valuelambda);
}
double GIS::getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY){
  assert(Feati<_sizeColValX && Featj<_sizeColValY);
  double lambda=_FM->getFeatureArraylambda(Feati,Featj, ilambdaX,ilambdaY);
  return lambda;
}
double**** GIS:: __getobs(){
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
double GIS::__getFeatconst(){
	double Featconst=1;

	int curr=0;
	for(int i=0; i< _sizeRowValX;i++){
		for(int j=0; j< _sizeY;j++){

			for(int delti=0; delti< _sizeColValX; delti++){
				for(int deltj=0; deltj< _sizeColValY; deltj++){
					for(int deltxi=0; deltxi < _sizeX; deltxi++){
						for(int deltyj=0; deltyj < _sizeY; deltyj++){
							for(int k=0; k< (*_FM).getMatrixIndexX(i,j).size();k++){
									curr++;
							}
						}
					}
					if(curr> Featconst) Featconst=curr;
					curr=0;
				}
			}
		}
	}
	return Featconst;
}
void GIS:: __getexp(){

  for(int i=0; i<_sizeColValX; i++){
    for( int j=0;j< _sizeColValY;j++){
      for(int k=0; k< _sizeX; k++){
        for(int l=0; l< _sizeY;l++){
          _expected[i][j][k][l]=0;
        }
      }
    }
  }
  for(int Feati=0; Feati< _sizeColValX; Feati++){
    for(int Featj=0; Featj< _sizeColValY; Featj++){

      for(int xi=0; xi< _sizeRowValX; xi++){
        _normaliser[Feati][Featj]=0;
        for(int yj=0; yj< _sizeY; yj++){
          _exponent[Feati][Featj][yj]=_FM->getFeatureArrayvalue(Feati,Featj,(*_valX)(xi,Feati), (*_Y)(yj,0));
          _normaliser[Feati][Featj]+=exp(_exponent[Feati][Featj][yj]);
        }
        for(int yj=0; yj< _sizeY; yj++){
          for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++){
            if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj){
              if((*_FM).getFeatureArraydelta(Feati,Featj,(*_FM).getMatrixIndexdX(xi,yj)[k],(*_FM).getMatrixIndexdY(xi,yj)[k],(*_valX)(xi,Feati),(*_Y)(yj,0))==1){
            	_expected[Feati][Featj][(*_FM).getMatrixIndexdX(xi,yj)[k]][(*_FM).getMatrixIndexdY(xi,yj)[k]]+=exp(_exponent[Feati][Featj][yj])/_normaliser[Feati][Featj];
              }
            }
          }
        }
      }
    }
  }
}
void GIS:: __gis(int maxit, double konv, bool test){
	/*
	 * 	time_t after3;
	time_t befor4;
	time_t after4;
	befor1=time(NULL);
	_gisTest=new GIS(*_valX,*_valY,*_X,*_Y,1,maxit,konv,false);
	after1=time(NULL);
	_timediff.push_back(difftime(after1,befor1));
	 */

	  //constant c for delta
	  double featconst = __getFeatconst();

	    for(int i=0; i<_sizeColValX; i++){
	      for(int j=0;j<_sizeColValY;j++){
	        for(int k=0;k< _sizeY;k++){
	          _exponent[i][j][k]=0;
	        }
	      }
	    }
	    for(int i=0; i< _sizeColValX; i++){
	      for(int j=0; j<_sizeColValY;j++){
	       _normaliser[i][j]=0;
	      }
	    }
	  int i=0;
	  double l=1;
	  double utime=0;
	  time_t befor;
	  time_t after;
	  while(utime<10 ){//&& fabs(l)>=konv
		befor=time(NULL);
	    l=0;
	    __getexp();
	    for(int Feati=0; Feati<_sizeColValX;Feati++){
	      for(int Featj=0; Featj< _sizeColValY;Featj++){
	        for(int lambdai=0; lambdai< _sizeX; lambdai++){
	          for(int lambdaj=0; lambdaj< _sizeY; lambdaj++){
	            double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,lambdai,lambdaj);
	            double newl=0;
	            if(fabs(_expected[Feati][Featj][lambdai][lambdaj]) < 0.00000001){_expected[Feati][Featj][lambdai][lambdaj]=0.01;}
	            if(fabs(_expected[Feati][Featj][lambdai][lambdaj]) > 0.00000001 && fabs(_observed[Feati][Featj][lambdai][lambdaj]) > 0.00000001 ){
	                    newl= oldl + (1/featconst)*log(_observed[Feati][Featj][lambdai][lambdaj]/_expected[Feati][Featj][lambdai][lambdaj]);
	            }
				else{
						newl=0;
				}
				(*_FM).setFeatureArraylambda(Feati,Featj,lambdai,lambdaj,newl);
				l+=fabs((_observed[Feati][Featj][lambdai][lambdaj]-_expected[Feati][Featj][lambdai][lambdaj]));
			   }
			 }
		   }
		 }
		 i++;
		 if(test){
			 _conv.push_back(l);
		 }
		  after=time(NULL);
		  utime+= difftime(after,befor);
	  	}
}
