#include "GISgp.h"

//training data, alphabete , startwert fuer lambda, startwert fuer delta, wert fuer sigma, test auf time, sekunden fuer den test
// GISgp::GISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,double lambdadeltaval, double sigma,int maxit,double konv, bool test,bool time,int seconds)
GISgp::GISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, IsParameter param)
  :IT(eX, eY, aX, aY, param, true)
{
  _param = param;
  _exponent = new double[_sizeY];
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
  if(param.time)
  {
    __gisgp(param.maxit, param.konv, param.lambdadeltaval, param.sigma, param.test, param.seconds);
  }
  else
  {
    __gisgp(param.maxit, param.konv, param.lambdadeltaval, param.sigma, param.test);
  }
}

GISgp::~GISgp()
{
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
  if(_delta != NULL){
    for(int i=0;i<_sizeColValX;i++){
      for(int j=0;j<_sizeColValY;j++){
        for(int k=0;k<_sizeX;k++){
          delete[] _delta[i][j][k];
        }
        delete[] _delta[i][j];
      }
      delete[] _delta[i];
    }
    delete[] _delta;
  }

  delete [] _exponent;

  _conv.clear();
}
void GISgp::__gisgp(int maxit, double konv, double lambdadeltaval, double sigma, bool test,int seconds){

  //constant c for delta
  double featconst = __getFeatconst();

  for(int k=0;k< _sizeY;k++){
    _exponent[k]=0;
  }
  _normaliser=0;

  // delta fuer die Newtonit.
  _delta= new double***[_sizeColValX];
  for(int i=0;i<_sizeColValX;i++)
  {
    _delta[i]= new double**[_sizeColValY];
    for(int j=0;j<_sizeColValY;j++)
    {
      _delta[i][j]= new double*[_sizeX];
      for(int k=0;k<_sizeX;k++)
      {
        _delta[i][j][k]= new double[_sizeY];
        for(int l=0;l<_sizeY;l++)
        {
          _delta[i][j][k][l]=lambdadeltaval;
        }
      }
    }
  }
  double l=1;
  double utime=0;
  time_t befor;
  time_t after;
  _iterations=0;
  while(utime<seconds ){
    befor=time(NULL);
    l=0;
    __getexp();
    for(int Feati=0; Feati<_sizeColValX;Feati++){
      for(int Featj=0; Featj< _sizeColValY;Featj++){
        for(int deltai=0; deltai< _sizeX; deltai++){
          for(int deltaj=0; deltaj< _sizeY; deltaj++){
            double newl;
            double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,deltai,deltaj);
            double z=1;
            while(fabs(z)>0.00001){
              z=(oldl+_delta[Feati][Featj][deltai][deltaj])/pow(sigma,2) +_expected[Feati][Featj][deltai][deltaj]*exp(_delta[Feati][Featj][deltai][deltaj]*featconst)- _observed[Feati][Featj][deltai][deltaj] ;
              double n= + 1/(pow(sigma,2))+_expected[Feati][Featj][deltai][deltaj]*featconst*exp(_delta[Feati][Featj][deltai][deltaj]*featconst);
              _delta[Feati][Featj][deltai][deltaj]= _delta[Feati][Featj][deltai][deltaj]-(z/n);
            }
            newl= oldl+_delta[Feati][Featj][deltai][deltaj];
            _FM->setFeatureArraylambda(Feati,Featj,deltai,deltaj,newl);
            l+=fabs((oldl+_delta[Feati][Featj][deltai][deltaj])/pow(sigma,2)+_expected[Feati][Featj][deltai][deltaj]*exp(_delta[Feati][Featj][deltai][deltaj]*featconst)- _observed[Feati][Featj][deltai][deltaj]);
          }
        }
      }
    }
    _iterations++;
    if(test){
      _conv.push_back(l);
    }
    after=time(NULL);
    utime+= difftime(after,befor);
  }
}
void GISgp::__gisgp(int maxit, double konv, double lambdadeltaval, double sigma, bool test){

  //constant c for delta
  double featconst = __getFeatconst();

  for(int i=0; i<_sizeColValX; i++){
    for(int j=0;j<_sizeColValY;j++){
      for(int k=0;k< _sizeY;k++){
        _exponent[k]=0;
      }
    }
  }
  // delta fuer die Newtonit.
  _delta= new double***[_sizeColValX];
  for(int i=0;i<_sizeColValX;i++){
    _delta[i]= new double**[_sizeColValY];
    for(int j=0;j<_sizeColValY;j++){
      _delta[i][j]= new double*[_sizeX];
      for(int k=0;k<_sizeX;k++){
        _delta[i][j][k]= new double[_sizeY];
        for(int l=0;l<_sizeY;l++){
          _delta[i][j][k][l]=lambdadeltaval;
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
            double newl;
            double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,lambdai,lambdaj);
            double z=1;
            while(z>0.00001){
              z=(oldl+_delta[Feati][Featj][lambdai][lambdaj])/pow(sigma,2) +_expected[Feati][Featj][lambdai][lambdaj]*exp(_delta[Feati][Featj][lambdai][lambdaj]*featconst)- _observed[Feati][Featj][lambdai][lambdaj] ;
              double n= + 1/(pow(sigma,2))+_expected[Feati][Featj][lambdai][lambdaj]*featconst*exp(_delta[Feati][Featj][lambdai][lambdaj]*featconst);
              _delta[Feati][Featj][lambdai][lambdaj]= _delta[Feati][Featj][lambdai][lambdaj]-(z/n);
            }
            newl= oldl+_delta[Feati][Featj][lambdai][lambdaj];
            _FM->setFeatureArraylambda(Feati,Featj,lambdai,lambdaj,newl);
            l+=fabs((oldl+_delta[Feati][Featj][lambdai][lambdaj])/pow(sigma,2)+_expected[Feati][Featj][lambdai][lambdaj]*exp(_delta[Feati][Featj][lambdai][lambdaj]*featconst)- _observed[Feati][Featj][lambdai][lambdaj]);
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

double GISgp:: getconv(int i){
  return _conv[i];
}
int GISgp:: getsizeconv(){
  return _conv.size();
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
        _normaliser=0;
        for(int yj=0; yj< _sizeY; yj++){
          _exponent[yj]=_FM->getFeatureArrayvalue(Feati,Featj,(*_valX)(xi,Feati), (*_Y)(yj,0));
          _normaliser+=exp(_exponent[yj]);
        }
        for(int yj=0; yj< _sizeY; yj++){
          for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++){
            if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj){
              if((*_FM).getFeatureArraydelta(Feati,Featj,(*_FM).getMatrixIndexdX(xi,yj)[k],(*_FM).getMatrixIndexdY(xi,yj)[k],(*_valX)(xi,Feati),(*_Y)(yj,0))==1){
                _expected[Feati][Featj][(*_FM).getMatrixIndexdX(xi,yj)[k]][(*_FM).getMatrixIndexdY(xi,yj)[k]]+=exp(_exponent[yj])/_normaliser;
              }
            }
          }
        }
      }
    }
  }
}
int GISgp:: getIterations(){
  return _iterations;
}
