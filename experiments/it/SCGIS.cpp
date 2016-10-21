#include "SCGIS.h"

#define EPSILON 0.00000001

// SCGIS::SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test,bool time,int seconds)
SCGIS::SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
:IT(eX, eY, aX, aY, systX, systY, param, false)
{
  _param = param;
  _exponent= new double**[_sizeSystX];
  for(int i=0;i<_sizeSystX; i++){
      _exponent[i]=new double*[_sizeRowValX];
      for(int xi=0;xi<_sizeRowValX;xi++){
        _exponent[i][xi]=new double[(int) pow(_Y->rows() ,_sizeColValY )];
        for(int y=0;y< pow(_Y->rows() ,_sizeColValY );y++){
          _exponent[i][xi][y]=0;
        }
      }

  }
  _normaliser=new double*[_sizeSystX];
  for(int i=0; i< _sizeSystX;i++){
      _normaliser[i]=new double[_sizeRowValX];
      for(int k=0;k<_sizeRowValX;k++){
        _normaliser[i][k]=pow(_Y->rows(),_sizeColValY);
      }
  }
  _delta=0.0;
  if(param.time)
  {
    __scgis(param.maxit,param.konv,param.test,param.seconds);
  }
  else
  {
    __scgis(param.maxit,param.konv,param.test);
  }
}

SCGIS::SCGIS(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
:IT(eX, eY, aX, aY, systX, systY, param, false)
{
  _param = param;
  _exponent= new double**[_sizeSystX];
  for(int i=0;i<_sizeSystX; i++){
      _exponent[i]=new double*[_sizeRowValX];
      for(int xi=0;xi<_sizeRowValX;xi++){
        _exponent[i][xi]=new double[(int) pow(_Y->rows() ,_sizeColValY )];
        for(int y=0;y< pow(_Y->rows() ,_sizeColValY );y++){
          _exponent[i][xi][y]=0;
        }
      }

  }
  _normaliser=new double*[_sizeSystX];
  for(int i=0; i< _sizeSystX;i++){
      _normaliser[i]=new double[_sizeRowValX];
      for(int k=0;k<_sizeRowValX;k++){
        _normaliser[i][k]=pow(_Y->rows(),_sizeColValY);
      }
  }
  _delta=0.0;
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
  for(int i=0;i<_sizeSystX;i++){
	  for(int j=0;j<_sizeRowValX;j++){
	      delete [] _exponent[i][j];
	  }
	  delete [] _exponent[i];
  }
  delete [] _exponent;

  for(int i=0;i<_sizeSystX;i++){
    delete [] _normaliser[i];
  }
  delete [] _normaliser;


  _conv.clear();

}


void SCGIS:: __scgis(int maxit, double konv, bool test,int seconds)
{
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
	  for(int feat=0;feat<_sizeSystX;feat++)
	  {
	    for(int delti=0;delti<pow(_sizeX,_systX[feat].size());delti++)
	      {
	        for(int deltj=0;deltj<pow(_sizeY,_systY[feat].size()); deltj++)
	        {
	          double expected = 0.0;
	          for(int y=0; y < pow(_sizeY,_sizeColValY); y++)
	          {
	            for(int k=0; k<_IM->getInstanceMatrixX(feat,delti,deltj).size();k++)
	            {
	              if(_IM->getInstanceMatrixY(feat,delti,deltj)[k] == y)
	              {
	                int x = _IM->getInstanceMatrixX(feat,delti,deltj)[k];
	                if(fabs(_normaliser[feat][x]) > EPSILON )
	                {
	                  // f_i(\bar{x}_j, y) is given by the 2 for and 1 if above
	                  expected+=exp(_exponent[feat][x][y])/_normaliser[feat][x];
	                }
	              }
	            }
	          }
	          double newl = 0.0;
	          if((fabs(expected)<EPSILON)  && (_observed[feat][delti][deltj]>EPSILON ))
	          {
	            expected=0.01;
	          }
	          if(fabs(_observed[feat][delti][deltj])>EPSILON)
	          {
	            _delta = log(_observed[feat][delti][deltj]/expected);
	            newl = _IM->getFeatureArraylambda(feat,delti,deltj) + _delta;
		        _IM->setFeatureArraylambda(feat,delti,deltj,newl);
	          }
	          else
	          {
	        	_delta = -1;               // evtl besseren Wert finden
	          }
	          l+=fabs(_observed[feat][delti][deltj]-expected);
	          for(int y=0;y< pow(_sizeY,_sizeColValY);y++)
	          {
	            for(int k=0; k<_IM->getInstanceMatrixX(feat,delti,deltj).size();k++)
	            {
	              if(_IM->getInstanceMatrixY(feat,delti,deltj)[k]==y)
	              {
	                int x=_IM->getInstanceMatrixX(feat,delti,deltj)[k];
	                _normaliser[feat][x]-=exp(_exponent[feat][x][y]);
	                _exponent[feat][x][y]+=_delta;
	                _normaliser[feat][x]+=exp(_exponent[feat][x][y]);
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
	//  cout << l << " " << _iterations << endl;
	  return l;
}
int SCGIS:: getIterations()
{
  return _iterations;
}

