#include "SCGIS.h"

#define EPSILON 0.00000001

// SCGIS::SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test,bool time,int seconds)
SCGIS::SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
:IT(eX, eY, aX, aY, systX, systY, param, false)
{
  _param = param;
  _exponent= new double**[_systX.size()];
  for(int i=0;i<_systX.size(); i++){
      _exponent[i]=new double*[_sizeRowValX];
      for(int xi=0;xi<_sizeRowValX;xi++){
        _exponent[i][xi]=new double[(int) pow(_Y->rows() ,_systY.size() )];
        for(int y=0;y<_sizeY;y++){
          _exponent[i][xi][y]=0;
        }
      }

  }
  _normaliser=new double*[_systX.size()];
  for(int i=0; i< _systX.size();i++){
      _normaliser[i]=new double[_sizeRowValX];
      for(int k=0;k<_sizeRowValX;k++){
        _normaliser[i][k]=pow(_Y->rows(),_systY.size());
      }
  }
  _delta= new double*[_sizeX];
  for(int i=0;i<_sizeX;i++){
    _delta[i]=new double[_sizeY];
    for(int j=0;j<_sizeY;j++){
      _delta[i][j]=0;
    }
  }
  if(param.time)
  {
    __scgis(param.maxit,param.konv,param.test,param.seconds);
  }
  else
  {
    __scgis(param.maxit,param.konv,param.test);
  }
}

double SCGIS:: getconv(int i){
  return _conv[i];
}

int SCGIS:: getsizeconv(){
  return _conv.size();
}

SCGIS::~SCGIS(){
  for(int i=0; i<_systX.size();i++){
      for(int k=0;k<(int) pow(_X->rows(),_sizeColValX);k++){
        delete [] _observed[i][k];
      }
    delete [] _observed[i];
  }
  delete [] _observed;


  for(int i=0;i<_systX.size();i++){
	  for(int j=0;j<_sizeRowValX;j++){
	      delete [] _exponent[i][j];
	  }
	  delete [] _exponent[i];
  }
  delete [] _exponent;

  for(int i=0;i<_systX.size();i++){
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
  for(int feat=0;feat<_systX.size();feat++)

      for(int delti=0;delti<pow(_X->rows(),_systX[feat].size());delti++)
      {
        for(int deltj=0;deltj<pow(_Y->rows(),_systY[feat].size()); deltj++)
        {
          double expected= 0.0;
          for(int y=0; y < pow(_Y->rows(),_sizeColValY); y++)
          {
            for(int k=0; k<_IM->getInstanceMatrixX(feat).size();k++)
            {
              if((_IM->getInstanceMatrixY(feat)[k] == y) && (_IM->getInstanceMatrixDeltaY(feat)[k] == deltj) && (_IM->getInstanceMatrixDeltaX(feat)[k] == delti))
              {
                int x = _IM->getInstanceMatrixX(feat)[k];
                cout << _normaliser[feat][x] << endl;
                if(fabs(_normaliser[feat][x]) > EPSILON )
                {
                   //f_i(\bar{x}_j, y) is given by the 2 for and 1 if above
                  expected+=exp(_exponent[feat][x][y])/_normaliser[feat][x];
                }
              }
            }
          }
          double newl = 0.0;
          if(fabs(expected)<EPSILON)
          {
            expected=0.01;
          }
          if(fabs(_observed[feat][delti][deltj])>EPSILON)
          {
            _delta[delti][deltj] = log(_observed[feat][delti][deltj]/expected);
            newl = _IM->getFeatureArraylambda(feat,delti,deltj) + _delta[delti][deltj];
          }
          else
          {
            newl = 0.0;
          }
          cout << _observed[feat][delti][deltj] << " " << expected << " " << newl  << " " << _delta[delti][deltj]<< endl;
          l+=fabs(_observed[feat][delti][deltj]-expected);
          _IM->setFeatureArraylambda(feat,delti,deltj,newl);
          for(int y=0;y<pow(_Y->rows(),_sizeColValY);y++)
          {
            for(int k=0; k<_IM->getInstanceMatrixX(feat).size();k++)
            {
              if((_IM->getInstanceMatrixY(feat)[k] == y) && (_IM->getInstanceMatrixDeltaY(feat)[k] == deltj) && (_IM->getInstanceMatrixDeltaX(feat)[k] == delti))
              {
                int x=_IM->getInstanceMatrixX(feat)[k];
                _normaliser[feat][x]-=exp(_exponent[feat][x][y]);
                _exponent[feat][x][y]+=_delta[delti][deltj];
                _normaliser[feat][x]+=exp(_exponent[feat][x][y]);
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
