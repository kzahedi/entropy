#include "SCGISgp.h"

#define EPSILON 0.00000001

SCGISgp::SCGISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
  :IT(eX, eY, aX, aY,systX,systY, param, false)

{
  _exponent= new double**[_sizeSystX];
  for(int i=0;i<_sizeSystX; i++)
  {
    _exponent[i]=new double*[_sizeRowValX];
    for(int xi=0;xi<_sizeRowValX;xi++)
    {
      _exponent[i][xi]=new double[(int) pow(_Y->rows(),_sizeColValY)];
      for(int y=0;y< pow(_Y->rows(),_sizeColValY);y++)
      {
        _exponent[i][xi][y]=0;
      }
    }
  }
  _normaliser=new double*[_sizeSystX];
  for(int i=0; i< _sizeSystX;i++)
  {
    _normaliser[i]=new double[_sizeRowValX];
    for(int k=0;k<_sizeRowValX;k++)
    {
      _normaliser[i][k]=pow(_Y->rows(),_sizeColValY);
    }
  }
  _delta= new double**[_sizeSystX];
  for(int i=0;i<_sizeSystX;i++)
  {
    _delta[i]= new double*[(int) pow(_X->rows(),_systX[i].size())];
    for(int k=0;k< pow(_X->rows(),_systX[i].size());k++)
    {
      _delta[i][k]= new double[(int) pow(_Y->rows(),_systY[i].size())];
      for(int l=0;l<pow(_Y->rows(),_systY[i].size());l++)
      {
        _delta[i][k][l]=param.lambdadeltaval;
      }
    }
  }
  if(param.time)
  {
    __scgis(param.maxit,param.konv,param.test,param.sigma,param.seconds);
  }
  else
  {
    __scgis(param.maxit,param.konv,param.test,param.sigma);
  }
}

SCGISgp:: ~SCGISgp()
{
  for(int i=0;i<_sizeSystX;i++)
  {
    for(int j=0;j< pow(_X->rows(),_sizeColValX);j++)
    {
      delete [] _exponent[i][j];
    }
    delete [] _exponent[i];
  }
  delete [] _exponent;

  for(int i=0;i<_sizeSystX;i++)
  {
    delete [] _normaliser[i];
  }
  delete [] _normaliser;

  for(int i=0;i<_sizeSystX;i++)
  {
    for(int j=0; j< pow(_X->rows(),_systX[i].size());j++)
    {
      delete [] _delta[i][j];
    }
    delete [] _delta[i];
  }
  delete [] _delta;

  _conv.clear();
}

double SCGISgp:: getconv(int i)
{
  return _conv[i];
}

int SCGISgp:: getsizeconv()
{
  return _conv.size();
}

double SCGISgp::__calculateIteration(bool test, double sigma)
{
  double l = 0.0;
  for(int feat=0; feat<_sizeSystX; feat++)
  {
    for(int delti=0; delti< pow(_X->rows(),_systX[feat].size());delti++)
	{
	  for(int deltj=0; deltj< pow(_Y->rows(), _systY[feat].size()); deltj++)
	  {
	    double expected=0;
	    for(int y=0;y < pow(_Y->rows(),_sizeColValY);y++)
	    {
	      for(int k=0; k<_IM->getInstanceMatrixX(feat,delti,deltj).size();k++)
	      {
	        if(_IM->getInstanceMatrixY(feat,delti,deltj)[k]==y)
	        {
	          int x=_IM->getInstanceMatrixX(feat,delti,deltj)[k];
	          if(fabs(_normaliser[feat][x])>0.00000001 )
	          {
	            expected+=exp(_exponent[feat][x][y])/_normaliser[feat][x];
	          }
	        }
	      }
	    }
	    double newl;
	    double oldl= (*_IM).getFeatureArraylambda(feat,delti,deltj);
        double zOld=2;
        double z=1;
        while(fabs(z-zOld)>0.0001)
        {
      	  zOld=z;
	      z=(oldl+_delta[feat][delti][deltj])/pow(sigma,2) +expected*exp(_delta[feat][delti][deltj])- _observed[feat][delti][deltj] ;
	      double n=  1/(pow(sigma,2))+expected*exp(_delta[feat][delti][deltj]);
	      _delta[feat][delti][deltj]= _delta[feat][delti][deltj]-(z/n);
	    }
	    newl= _IM->getFeatureArraylambda(feat,delti,deltj)+_delta[feat][delti][deltj];
	    l+=fabs((oldl+_delta[feat][delti][deltj])/pow(sigma,2) +expected*exp(_delta[feat][delti][deltj])- _observed[feat][delti][deltj] );
	    _IM->setFeatureArraylambda(feat,delti,deltj,newl);
	    for(int y=0;y<pow(_Y->rows(),_sizeColValY);y++)
	    {
	      for(int k=0; k<_IM->getInstanceMatrixX(feat,delti,deltj).size();k++)
	      {
	        if(_IM->getInstanceMatrixY(feat,delti,deltj)[k]==y)
	        {
	          int x=_IM->getInstanceMatrixX(feat,delti,deltj)[k];
	          _normaliser[feat][x]-=exp(_exponent[feat][x][y]);
	          _exponent[feat][x][y]+=_delta[feat][delti][deltj];
	          _normaliser[feat][x]+=exp(_exponent[feat][x][y]);
	        }
	      }
	    }
	  }
	}
  }
  _iterations++;
  if(test)
  {
    _conv.push_back(l);
  }
  return l;
}
void SCGISgp:: __scgis(int maxit, double konv,bool test,double sigma,int seconds)
{
  double utime=0;
  _iterations=0;
  time_t befor;
  time_t after;
  while(utime<seconds)
  {
    befor=time(NULL);
	__calculateIteration(test,sigma);
	after=time(NULL);
	utime+= difftime(after,befor);
  }
}

void SCGISgp:: __scgis(int maxit, double konv,bool test,double sigma)
{
	cout << " scgis " << endl;
  double l=1;
  while(_iterations<maxit  ) //&& fabs(l)>konv
  {
    l=__calculateIteration(test,sigma);
    cout << l << endl;
  }
}

int SCGISgp:: getIterations()
{
  return _iterations;
}
