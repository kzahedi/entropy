#include "GIS.h"

GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv) {
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
      _FM=new FeatureMatrix(*_valX,*_valY,*_X,*_Y,lambdavalue);
      __gis(maxit, konv);
}
// ColValY - Anzahl der Y-Knoten
GIS::GIS(int ColValY, DContainer &eX){
    _sizeX=2;
    _sizeY=2;
    _valX= &eX;
    _sizeColValY= ColValY;
    _sizeColValX= (*_valX).columns();
    _sizeRowValX= (*_valX).rows();
    _valY= new DContainer(_sizeRowValX,ColValY);
    _sizeRowValY= 0;
    _X=new DContainer(_sizeX,1);
    (*_X) << 0 << 1;
    _Y=new DContainer(_sizeY,1);
    (*_Y) << 0 << 1;
    _FM=new FeatureMatrix(*_valX,*_valY,*_X,*_Y,0);
}
double GIS::gis(int Feati,int Featj,double ValX,double ValY){
  double norm=0;
  double exponent=0;
    exponent= exp((*_FM).getFeatureArrayvalueval(Feati,Featj,ValX,ValY) );
  for(int yi=0;yi<_sizeY;yi++){
    norm+= exp( (*_FM).getFeatureArrayvalueval(Feati,Featj,ValX,(*_Y)(yi,0)));
  }
  return exponent/norm;
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
  double**** observed;
    observed = new double***[_sizeColValX];
    for(int i=0; i<_sizeColValX; i++){
      observed[i]=new double**[_sizeColValY];
      for( int j=0;j< _sizeColValY;j++){
        observed[i][j]=new double*[_sizeX];
        for(int k=0; k< _sizeX; k++){
          observed[i][j][k]= new double[_sizeY];
          for(int l=0; l< _sizeY;l++){
            observed[i][j][k][l]=0;
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
            	if(_FM->getFeatureArraydeltaval(Feati,Featj,delti,deltj,(*_valX)(i,Feati),(*_valY)(i,Featj))!=-1){
            		observed[Feati][Featj][delti][deltj]++;
            	}
            }
          }
        }
      }
    }
    return observed;


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
void GIS:: __getexp(double**** &expect, double*** &exponent,double** &normaliser){

  for(int i=0; i<_sizeColValX; i++){
    for( int j=0;j< _sizeColValY;j++){
      for(int k=0; k< _sizeX; k++){
        for(int l=0; l< _sizeY;l++){
          expect[i][j][k][l]=0;
        }
      }
    }
  }
  for(int Feati=0; Feati< _sizeColValX; Feati++){
    for(int Featj=0; Featj< _sizeColValY; Featj++){

      for(int xi=0; xi< _sizeRowValX; xi++){
        normaliser[Feati][Featj]=0;
        for(int yj=0; yj< _sizeY; yj++){
          exponent[Feati][Featj][yj]=_FM->getFeatureArrayvalueval(Feati,Featj,(*_valX)(xi,Feati), (*_Y)(yj,0));
          normaliser[Feati][Featj]+=exp(exponent[Feati][Featj][yj]);
       // cout << exp(exponent[Feati][Featj][yj]) << endl;
        }
        for(int yj=0; yj< _sizeY; yj++){
          for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++){
            if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj){
            	if((*_FM).getFeatureArraydelta(Feati,Featj,(*_FM).getMatrixIndexdX(xi,yj)[k],(*_FM).getMatrixIndexdY(xi,yj)[k],xi,yj)==1){
            		expect[Feati][Featj][(*_FM).getMatrixIndexdX(xi,yj)[k]][(*_FM).getMatrixIndexdY(xi,yj)[k]]+=exp(exponent[Feati][Featj][yj])/normaliser[Feati][Featj];
            	}
            }
          }
        }
      }
    }
  }
}
void GIS:: __normaliselambdafeat(){
  double** sumlambda;
  sumlambda=new double*[_sizeColValX];
  for(int i=0;i<_sizeColValX;i++){
    sumlambda[i]=new double[_sizeColValX];
    for(int j=0;j<_sizeColValY;j++){
      sumlambda[i][j]=0;
    }
  }
  double newl=0;
  for(int Feati=0;Feati< _sizeColValX;Feati++){
    for(int Featj=0;Featj< _sizeColValY;Featj++){
      for(int deltai=0;deltai< _sizeX; deltai++){
        for(int deltaj=0;deltaj<_sizeY;deltaj++){
          sumlambda[Feati][Featj]+=fabs(_FM->getFeatureArraylambda(Feati,Featj,deltai,deltaj));
        }
      }
    }
  }
  for(int Feati=0;Feati< _sizeColValX;Feati++){
    for(int Featj=0;Featj< _sizeColValY;Featj++){
      for(int deltai=0;deltai< _sizeX; deltai++){
        for(int deltaj=0;deltaj<_sizeY;deltaj++){
          newl=_FM->getFeatureArraylambda(Feati,Featj,deltai,deltaj)/sumlambda[Feati][Featj];
          _FM->setFeatureArraylambda(Feati,Featj,deltai,deltaj,newl);
        }
      }
    }
  }

}
void GIS:: __normaliselambda(){
  double sumlambda;
  double newl=0;
  for(int Feati=0;Feati< _sizeColValX;Feati++){
    for(int Featj=0;Featj< _sizeColValY;Featj++){
      for(int deltai=0;deltai< _sizeX; deltai++){
        for(int deltaj=0;deltaj<_sizeY;deltaj++){
          sumlambda+=fabs(_FM->getFeatureArraylambda(Feati,Featj,deltai,deltaj));
        }
      }
    }
  }
  for(int Feati=0;Feati< _sizeColValX;Feati++){
    for(int Featj=0;Featj< _sizeColValY;Featj++){
      for(int deltai=0;deltai< _sizeX; deltai++){
        for(int deltaj=0;deltaj<_sizeY;deltaj++){
          newl=_FM->getFeatureArraylambda(Feati,Featj,deltai,deltaj)/sumlambda;
          _FM->setFeatureArraylambda(Feati,Featj,deltai,deltaj,newl);
        }
      }
    }
  }

}
void GIS:: __gis(int maxit, double konv){
  //observed
  double**** observ = __getobs();

  //constant c for delta
  double featconst = __getFeatconst();

  //array for exp
  double**** expected;
    expected = new double***[_sizeColValX];
    for(int i=0; i<_sizeColValX; i++){
      expected[i]=new double**[_sizeColValY];
      for( int j=0;j< _sizeColValY;j++){
        expected[i][j]=new double*[_sizeX];
        for(int k=0; k< _sizeX; k++){
          expected[i][j][k]= new double[_sizeY];
          for(int l=0; l< _sizeY;l++){
            expected[i][j][k][l]=0;
          }
        }
      }
    }
  double*** exponent;
    exponent = new double**[_sizeColValX];
    for(int i=0; i<_sizeColValX; i++){
      exponent[i]=new double*[_sizeColValY];
      for(int j=0;j<_sizeColValY;j++){
        exponent[i][j]=new double[_sizeRowValY];
        for(int k=0;k< _sizeY;k++){
          exponent[i][j][k]=0;
        }
      }
    }
  double** normaliser;
    normaliser= new double*[_sizeColValX];
    for(int i=0; i< _sizeColValX; i++){
      normaliser[i]=new double[_sizeColValY];
      for(int j=0; j<_sizeColValY;j++){
        normaliser[i][j]=0;
      }
    }
  int i=0;
  double l=1;
  bool norm=false;
  while(i<maxit && fabs(l)>=konv ){
    l=0;
    __getexp(expected,exponent,normaliser);
    for(int Feati=0; Feati<_sizeColValX;Feati++){
      for(int Featj=0; Featj< _sizeColValY;Featj++){
        for(int lambdai=0; lambdai< _sizeX; lambdai++){
          for(int lambdaj=0; lambdaj< _sizeY; lambdaj++){
            double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,lambdai,lambdaj);
            double newl=0;
            //cout << "exp" << expected[Feati][Featj][lambdai][lambdaj] << endl;
            if(fabs(expected[Feati][Featj][lambdai][lambdaj]) < 0.00000001){ expected[Feati][Featj][lambdai][lambdaj]=0.01;}
            //cout << "exp " <<  expected[Feati][Featj][lambdai][lambdaj] << endl;
            //cout << observ[Feati][Featj][lambdai][lambdaj] << endl;
            if(fabs(expected[Feati][Featj][lambdai][lambdaj]) > 0.0000001 && fabs(observ[Feati][Featj][lambdai][lambdaj]) > 0.00000001 ){
                    double p=(observ[Feati][Featj][lambdai][lambdaj]/expected[Feati][Featj][lambdai][lambdaj]);
                    //cout << "obs " <<  observ[Feati][Featj][lambdai][lambdaj] << endl;
                    //cout << expected[Feati][Featj][lambdai][lambdaj] << endl;
                    //cout << p << endl;
                    //cout << log(p) << endl;
                   // if (lambdai == 0 && lambdaj == 0)
                   // {
                    //  cout << observ[Feati][Featj][lambdai][lambdaj] << " " <<
                     //  expected[Feati][Featj][lambdai][lambdaj] << endl;
                     // cout << "o: " << oldl << " fc: " << featconst << " p: " << p << " l: " << log(p);
                   // }
                    newl= oldl + (1/featconst)*log(p);
                  //  if (lambdai == 0 && lambdaj == 0)
                    //  cout << " n: " << newl << endl;
                    //cout << "oldl "<< oldl << " "<< newl << endl;

		 							}
						else{newl=0;}
							(*_FM).setFeatureArraylambda(Feati,Featj,lambdai,lambdaj,newl);
							l+=fabs((observ[Feati][Featj][lambdai][lambdaj]-expected[Feati][Featj][lambdai][lambdaj]));
							if(fabs(newl>100)){
							//__normaliselambdafeat();
							//__normaliselambda();
							norm =true;
							}
					}
				}
			}
		}
		i++;
		cout << l << endl;
	}
	cout << norm << endl;
}

