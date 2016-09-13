#include "GISgp.h"

//training data, alphabete , startwert fuer lambda, startwert fuer delta, wert fuer sigma
GISgp::GISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,double lambdadeltaval, double sigma,int maxit,double konv, bool test) {
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
	  __gisgp(maxit, konv ,lambdadeltaval,sigma,test);
}
void GISgp::__gisgp(int maxit, double konv, double lambdadeltaval, double sigma, bool test){

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
	   // delta fuer die Newtonit.
	  _lambdadelta= new double***[_sizeColValX];
	  for(int i=0;i<_sizeColValX;i++){
		  _lambdadelta[i]= new double**[_sizeColValY];
		  for(int j=0;j<_sizeColValY;j++){
			  _lambdadelta[i][j]= new double*[_sizeX];
			  for(int k=0;k<_sizeX;k++){
				  _lambdadelta[i][j][k]= new double[_sizeY];
				  for(int l=0;l<_sizeY;l++){
					  _lambdadelta[i][j][k][l]=lambdadeltaval;
				  }
			  }
		  }
	  }
	  int i=0;
	  double l=1;
	  while(i<maxit && fabs(l)>=konv ){
	    l=0;
	    __getexp();
	    for(int Feati=0; Feati<_sizeColValX;Feati++){
	      for(int Featj=0; Featj< _sizeColValY;Featj++){
	        for(int lambdai=0; lambdai< _sizeX; lambdai++){
	          for(int lambdaj=0; lambdaj< _sizeY; lambdaj++){
	        	//_lambdadelta[Feati][Featj][lambdai][lambdaj]=lambdadeltaval;
	        	double newl;
	            double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,lambdai,lambdaj);
	            double z=1;
	          //  cout << "observed 1 " << _observed[Feati][Featj][lambdai][lambdaj] << "  expected " << (oldl+lambdadelta[Feati][Featj][lambdai][lambdaj])/pow(sigma,2)+_expected[Feati][Featj][lambdai][lambdaj]*exp(lambdadelta[Feati][Featj][lambdai][lambdaj]*featconst) << endl;
	               while(z>0.001){
	            	   z=(oldl+_lambdadelta[Feati][Featj][lambdai][lambdaj])/pow(sigma,2) +_expected[Feati][Featj][lambdai][lambdaj]*exp(_lambdadelta[Feati][Featj][lambdai][lambdaj]*featconst)- _observed[Feati][Featj][lambdai][lambdaj] ;
	            	   double n= + 1/(pow(sigma,2))+_expected[Feati][Featj][lambdai][lambdaj]*featconst*exp(_lambdadelta[Feati][Featj][lambdai][lambdaj]*featconst);
	            	   _lambdadelta[Feati][Featj][lambdai][lambdaj]= _lambdadelta[Feati][Featj][lambdai][lambdaj]-(z/n);
	            	//   cout <<lambdadelta[Feati][Featj][lambdai][lambdaj] << endl;
	         //   	   cout << (oldl+lambdadelta[Feati][Featj][lambdai][lambdaj])/pow(sigma,2)- _observed[Feati][Featj][lambdai][lambdaj]+ _expected[Feati][Featj][lambdai][lambdaj]*exp(lambdadelta[Feati][Featj][lambdai][lambdaj]*featconst)<< endl;
	               }
	             //  cout << lambdadelta[Feati][Featj][lambdai][lambdaj] << endl;
	               newl= oldl+_lambdadelta[Feati][Featj][lambdai][lambdaj];
	              // cout <<"newl " << newl << " " << oldl << endl;
	            _FM->setFeatureArraylambda(Feati,Featj,lambdai,lambdaj,newl);
	            l+=fabs((oldl+_lambdadelta[Feati][Featj][lambdai][lambdaj])/pow(sigma,2)+_expected[Feati][Featj][lambdai][lambdaj]*exp(_lambdadelta[Feati][Featj][lambdai][lambdaj]*featconst)- _observed[Feati][Featj][lambdai][lambdaj]);
	       //     cout << (oldl+lambdadelta[Feati][Featj][lambdai][lambdaj])/pow(sigma,2)- _observed[Feati][Featj][lambdai][lambdaj]+ _expected[Feati][Featj][lambdai][lambdaj]*exp(lambdadelta[Feati][Featj][lambdai][lambdaj]*featconst)<< endl;

	            //cout << "observed " << _observed[Feati][Featj][lambdai][lambdaj]-(oldl+lambdadelta[Feati][Featj][lambdai][lambdaj])/pow(sigma,2) << "  expected " << _expected[Feati][Featj][lambdai][lambdaj] << endl;

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
double GISgp::prop(int Feati,int Featj,double ValX,double ValY){
  double norm=0;
  double exponent=0;
  exponent= exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,ValY) );
  for(int yi=0;yi<_sizeY;yi++){
	    norm+= exp( (*_FM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
  }
  return exponent/norm;
}
double GISgp::prop(int rowX,vector<vector<double> > Y, int rowY){
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
double GISgp:: prop(vector<vector<double> > X,int rowX,vector<vector<double> > Y, int rowY){
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
double GISgp:: getconv(int i){
	return _conv[i];
}
int GISgp:: getsizeconv(){
	return _conv.size();
}
void GISgp::setFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY,double valuelambda){
  assert(Feati<_sizeColValX && Featj<_sizeColValY);
  _FM->setFeatureArraylambda(Feati,Featj,ilambdaX,ilambdaY,valuelambda);
}
double GISgp::getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY){
  assert(Feati<_sizeColValX && Featj<_sizeColValY);
  double lambda=_FM->getFeatureArraylambda(Feati,Featj, ilambdaX,ilambdaY);
  return lambda;
}
double**** GISgp:: __getobs(){
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
double GISgp::__getFeatconst(){
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
void GISgp:: __getexp(){

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
