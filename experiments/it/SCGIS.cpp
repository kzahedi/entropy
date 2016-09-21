#include "SCGIS.h"

#define EPSILON 0.00000001

SCGIS::SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test,bool time,int seconds)
  :IT(eX, eY, aX, aY, lambdavalue,false)
{

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

  _delta= new double*[_sizeX];
  for(int i=0;i<_sizeX;i++){
    _delta[i]=new double[_sizeY];
    for(int j=0;j<_sizeY;j++){
      _delta[i][j]=0;
    }
  }

  if(time)
  {
    __scgis(maxit,konv,test,seconds);
  }
  else
  {
    __scgis(maxit,konv,test);
  }
}

double SCGIS:: getconv(int i){
  return _conv[i];
}

int SCGIS:: getsizeconv(){
  return _conv.size();
}

SCGIS::~SCGIS(){

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


void SCGIS:: __scgis(int maxit, double konv, bool test,int seconds){
  double l=1;
  double utime=0;
  time_t befor;
  time_t after;
  _iterations = 0;
  while(utime<seconds ){
    befor=time(NULL);
    __calculateIteration(test);
    after=time(NULL);
    utime+= difftime(after,befor);
  }

}

void SCGIS:: __scgis(int maxit, double konv, bool test)
{
  double l=1;
  _iterations = 0;
  while(_iterations < maxit && fabs(l)>konv ) 
  {
    l = __calculateIteration(test);
  }
}

double SCGIS::__calculateIteration(bool test)
{
  double l=0;
  for(int Feati=0;Feati<_sizeColValX;Feati++)
  {
    for(int Featj=0;Featj<_sizeColValY;Featj++)
    {
      for(int delti=0;delti<_sizeX;delti++)
      {
        for(int deltj=0;deltj<_sizeY; deltj++)
        {
          double expected = 0.0;
          for(int y=0; y < _sizeY; y++)
          {
            for(int k=0; k<_IM->getInstanceMatrixX(Feati,Featj,delti,deltj).size();k++)
            {
              if(_IM->getInstanceMatrixY(Feati,Featj,delti,deltj)[k] == y)
              {
                int x = _IM->getInstanceMatrixX(Feati,Featj,delti,deltj)[k];
                if(fabs(_normaliser[Feati][Featj][x]) > EPSILON )
                {
                  // f_i(\bar{x}_j, y) is given by the 2 for and 1 if above
                  expected+=exp(_exponent[Feati][Featj][x][y])/_normaliser[Feati][Featj][x];
                }
              }
            }
          }
          double newl = 0.0;
          if(fabs(expected)<EPSILON)
          {
            expected=0.01;
          }
          if(fabs(_observed[Feati][Featj][delti][deltj])>EPSILON)
          {
            _delta[delti][deltj] = log(_observed[Feati][Featj][delti][deltj]/expected);
            newl = _IM->getFeatureArraylambda(Feati,Featj,delti,deltj) + _delta[delti][deltj];
          }
          else
          {
            newl = 0.0;
          }
          l+=fabs(_observed[Feati][Featj][delti][deltj]-expected);
          _IM->setFeatureArraylambda(Feati,Featj,delti,deltj,newl);
          for(int y=0;y<_sizeY;y++)
          {
            for(int k=0; k<_IM->getInstanceMatrixX(Feati,Featj,delti,deltj).size();k++)
            {
              if(_IM->getInstanceMatrixY(Feati,Featj,delti,deltj)[k]==y)
              {
                int x=_IM->getInstanceMatrixX(Feati,Featj,delti,deltj)[k];
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
  _iterations++;
  if(test){
    _conv.push_back(l);
  }
  return l;
}

int SCGIS:: getIterations(){
  return _iterations;
}
