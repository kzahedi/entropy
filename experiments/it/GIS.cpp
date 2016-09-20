#include "GIS.h"

GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test, bool time, int seconds)
:IT(eX, eY, aX, aY, lambdavalue,true)
{
  _exponent   = new double[_sizeRowValY];
  _normaliser = 0.0;
  _expected   = new double***[_sizeColValX];
  for(int i=0; i<_sizeColValX; i++)
  {
    _expected[i]=new double**[_sizeColValY];
    for( int j=0;j< _sizeColValY;j++)
    {
      _expected[i][j]=new double*[_sizeX];
      for(int k=0; k< _sizeX; k++)
      {
        _expected[i][j][k]= new double[_sizeY];
        for(int l=0; l< _sizeY;l++)
        {
          _expected[i][j][k][l]=0;
        }
      }
    }
  }

  if(time)
  {
    __gis(maxit, konv, test, seconds);
  }
  else
  {
    __gis(maxit, konv, test);
  }
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
  if(_exponent != NULL){
    delete [] _exponent;
  }
  _conv.clear();
}
double GIS:: getconv(int i){
  return _conv[i];
}
int GIS:: getsizeconv(){
  return _conv.size();
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
  for(int Feati=0; Feati< _sizeColValX; Feati++)
  {
    for(int Featj=0; Featj< _sizeColValY; Featj++)
    {
      for(int xi=0; xi< _sizeRowValX; xi++)
      {
        _normaliser = 0.0;
        for(int yj=0; yj< _sizeY; yj++)
        {
          _exponent[yj]=_FM->getFeatureArrayvalue(Feati,Featj,(*_valX)(xi,Feati), (*_Y)(yj,0));
          _normaliser += exp(_exponent[yj]);
        }
        for(int yj=0; yj< _sizeY; yj++)
        {
          for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++)
          {
            if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj)
            {
              if((*_FM).getFeatureArraydelta(Feati,
                                             Featj,
                                             (*_FM).getMatrixIndexdX(xi,yj)[k],
                                             (*_FM).getMatrixIndexdY(xi,yj)[k],
                                             (*_valX)(xi,Feati),
                                             (*_Y)(yj,0)) == 1)
              {
              _expected[Feati]
                       [Featj]
                       [(*_FM).getMatrixIndexdX(xi,yj)[k]]
                       [(*_FM).getMatrixIndexdY(xi,yj)[k]] += exp(_exponent[yj])/_normaliser;
              }
            }
          }
        }
      }
    }
  }
}
void GIS:: __gis(int maxit, double konv, bool test,int seconds){

    //constant c for delta
    double featconst = __getFeatconst();

    for(int k=0;k< _sizeY;k++){
      _exponent[k]=0;
    }
    _normaliser=0;
    double l=1;
    double utime=0;
    time_t befor;
    time_t after;
    _iterations=0;
    while(utime<seconds )
    {
      befor=time(NULL);
      l=0;
      __getexp();
      for(int Feati=0; Feati<_sizeColValX;Feati++){
        for(int Featj=0; Featj< _sizeColValY;Featj++){
          for(int deltai=0; deltai< _sizeX; deltai++){
            for(int deltaj=0; deltaj< _sizeY; deltaj++){
              double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,deltai,deltaj);
              double newl=0;

              if(fabs(_expected[Feati][Featj][deltai][deltaj]) < 0.00000001)
              {
                _expected[Feati][Featj][deltai][deltaj]=0.01;
              }
              if(fabs(_expected[Feati][Featj][deltai][deltaj]) > 0.00000001 &&
                 fabs(_observed[Feati][Featj][deltai][deltaj]) > 0.00000001 )
              {
                newl = oldl + (1/featconst)*log(_observed[Feati][Featj][deltai][deltaj]/_expected[Feati][Featj][deltai][deltaj]);
              }
              else
              {
                newl=0;
              }
              (*_FM).setFeatureArraylambda(Feati,Featj,deltai,deltaj,newl);
              l += fabs(_observed[Feati][Featj][deltai][deltaj] - _expected[Feati][Featj][deltai][deltaj]);
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
void GIS:: __gis(int maxit, double konv, bool test){

    //constant c for delta
    double featconst = __getFeatconst();

    for(int k=0;k< _sizeY;k++){
      _exponent[k]=0;
    }
    _normaliser=0;
    int i=0;
    double l=1;
    while(i<maxit && fabs(l)>=konv ){
      l=0;
      __getexp();
      for(int Feati=0; Feati<_sizeColValX;Feati++){
        for(int Featj=0; Featj< _sizeColValY;Featj++){
          for(int deltai=0; deltai< _sizeX; deltai++){
            for(int deltaj=0; deltaj< _sizeY; deltaj++){
              double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,deltai,deltaj);
              double newl=0;
              if(fabs(_expected[Feati][Featj][deltai][deltaj]) < 0.00000001){_expected[Feati][Featj][deltai][deltaj]=0.01;}
              if(fabs(_expected[Feati][Featj][deltai][deltaj]) > 0.00000001 && fabs(_observed[Feati][Featj][deltai][deltaj]) > 0.00000001 ){
                      newl= oldl + (1/featconst)*log(_observed[Feati][Featj][deltai][deltaj]/_expected[Feati][Featj][deltai][deltaj]);
              }
        else{
            newl=0;
        }
        (*_FM).setFeatureArraylambda(Feati,Featj,deltai,deltaj,newl);
        l+=fabs((_observed[Feati][Featj][deltai][deltaj]-_expected[Feati][Featj][deltai][deltaj]));
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
int GIS:: getIterations(){
  return _iterations;
}

