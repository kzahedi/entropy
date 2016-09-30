#include "GISgp.h"

//training data, alphabete , startwert fuer lambda, startwert fuer delta, wert fuer sigma, test auf time, sekunden fuer den test
// GISgp::GISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,double lambdadeltaval, double sigma,int maxit,double konv, bool test,bool time,int seconds)
GISgp::GISgp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
  :IT(eX, eY, aX, aY, systX, systY, param, true)
{
  _normaliser = new double[_sizeSystX];

  _exponent= new double*[_sizeSystX];
  for(int i=0;i<_sizeSystX;i++){
	  _exponent[i]= new double[(int) pow(_Y->rows(),_sizeColValY)];
  }
  _expected = new double**[_sizeSystX];
  for(int i=0; i<_sizeSystX; i++)
  {
    _expected[i]=new double*[(int) pow(_X->rows(),_systX[i].size())];
    for(int k=0; k< (int) pow(_X->rows(),_systX[i].size()); k++)
    {
      _expected[i][k]= new double[(int) pow(_Y->rows(),_systY[i].size())];
      for(int l=0; l< _sizeY;l++)
      {
        _expected[i][k][l]=0;
      }
    }
  }
  _delta= new double**[_sizeSystX];
  for(int i=0;i<_sizeSystX;i++)
  {
    _delta[i]= new double*[(int) pow(_X->rows(),_systX[i].size())];
    for(int j=0;j<(int) pow(_X->rows(),_systX[i].size());j++)
    {
      _delta[i][j]= new double[(int) pow(_Y->rows(),_systY[i].size())];
      for(int k=0;k<(int) pow(_Y->rows(),_systY[i].size());k++)
      {
        _delta[i][j][k]=param.lambdadeltaval;
      }
    }
  }
  if(param.time)
  {
    __gisgp(param.maxit, param.konv, param.sigma, param.test, param.seconds);
  }
  else
  {
    __gisgp(param.maxit, param.konv, param.sigma, param.test);
  }
}

GISgp::~GISgp()
{
  if(_expected != NULL)
  {
    for(int i=0;i<_sizeSystX;i++)
    {
      for(int k=0;k<pow(_X->rows(),_systX[i].size());k++)
      {
        delete [] _expected[i][k];
      }
      delete[] _expected[i];
    }
    delete[] _expected;
  }

  if(_delta != NULL)
  {
    for(int i=0;i<_sizeSystX;i++)
    {
      for(int j=0;j<pow(_X->rows(),_systX[i].size());j++)
      {
        delete[] _delta[i][j];
      }
      delete[] _delta[i];
    }
    delete[] _delta;
  }

  delete [] _exponent;

  _conv.clear();
  delete [] _normaliser;
}
double GISgp::__calculateIteration(double featconst,double sigma, bool test)
{
	    double l=0;
	    __getexp();
	    for(int feat=0; feat<_sizeSystX; feat++)
	    {
          for(int deltai=0; deltai<  pow(_X->rows(),_systX[feat].size()); deltai++)
	      {
	        for(int deltaj=0; deltaj<  pow(_Y->rows(),_systY[feat].size()); deltaj++)
	        {
	          double newl;
	          double oldl= (*_FM).getFeatureArraylambda(feat,deltai,deltaj);
	          double zOld=2;
	          double z=1;
	          while(fabs(z-zOld)>0.0001)
	          {
	        	zOld=z;
	            z=(oldl+_delta[feat][deltai][deltaj])/pow(sigma,2) +_expected[feat][deltai][deltaj]*exp(_delta[feat][deltai][deltaj]*featconst)- _observed[feat][deltai][deltaj] ;
	            double n= + 1/(pow(sigma,2))+_expected[feat][deltai][deltaj]*featconst*exp(_delta[feat][deltai][deltaj]*featconst);
	            if(fabs(n)<0.00000001){ cout << "  hieer " << endl;}
	            _delta[feat][deltai][deltaj]= _delta[feat][deltai][deltaj]-(z/n);
	          }
	          newl= oldl+_delta[feat][deltai][deltaj];
	          _FM->setFeatureArraylambda(feat,deltai,deltaj,newl);
	          l+=fabs((oldl+_delta[feat][deltai][deltaj])/pow(sigma,2)+_expected[feat][deltai][deltaj]*exp(_delta[feat][deltai][deltaj]*featconst)- _observed[feat][deltai][deltaj]);
	         }
	       }
	     }
	    _iterations++;
	    if(test){
	      _conv.push_back(l);
	    }
	 //   cout << "l" <<l << " " << _iterations << endl;
        return l;
}
void GISgp::__gisgp(int maxit, double konv, double sigma, bool test,int seconds)
{
  //constant c for delta
  double featconst = __getFeatconst();

  for(int i=0;i<_sizeSystX;i++)
  {
	  for(int j=0;j<pow(_Y->rows(),_sizeColValY);j++)
	  {
		  _exponent[i][j]= 0.0;
	  }
  }
  double l=1;
  double utime=0;
  time_t befor;
  time_t after;
  _iterations=0;
  while(utime<seconds)
  {
  befor=time(NULL);
  l=0;
   __calculateIteration(featconst,sigma,test);
   after=time(NULL);
   utime+= difftime(after,befor);
  }
}

void GISgp::__gisgp(int maxit, double konv, double sigma, bool test)
{
  //constant c for delta
  double featconst = __getFeatconst();

  for(int i=0;i<_sizeSystX;i++)
  {
	  for(int j=0;j<pow(_Y->rows(),_sizeColValY);j++)
	  {
		  _exponent[i][j]= 0.0;
	  }
  }

  double l=1;
  _iterations=0;
  while(_iterations<maxit && fabs(l)>=konv )
  {
	  l=__calculateIteration(featconst,sigma,true);
  }
}

double GISgp:: getconv(int i)
{
  return _conv[i];
}

int GISgp:: getsizeconv()
{
  return _conv.size();
}

double GISgp::__getFeatconst()
{
  double r = 0.0;
  for(int i=0; i< _sizeRowValX;i++) // i-th data row
  {
    for(int j=0; j< pow(_Y->rows(),_sizeColValY);j++) // y-alphabet
	{
	  int v = (*_FM).getMatrixIndexFeat(i,j).size(); // the number of matching deltas
	  if(v > r) r = v;
	}
  }
  return r;
}

void GISgp::__getexp()
{
	  for(int i=0; i<_sizeSystX; i++)
	  {
	    for(int k=0; k< (int) pow(_X->rows(),_systX[i].size()); k++)
	    {
	      for(int l=0; l< (int) pow(_Y->rows(),_systY[i].size());l++)
	      {
	        _expected[i][k][l]=0;
	      }
	    }
	  }

	  for(int xi=0; xi< _sizeRowValX; xi++)
	  {
	    for(int i=0;i<_sizeSystX;i++){
	      _normaliser[i]=0;
	    }
	    for(int yj=0; yj< pow(_Y->rows(),_sizeColValY); yj++)
	    {
	      for(int k=0; k< (*_FM).getMatrixIndexFeat(xi,yj).size();k++){
	    	int index=_FM->getMatrixIndexFeat(xi,yj)[k];
	    	_exponent[index][yj]=_FM->getFeatureArrayvalueAlphY(index,xi,yj);
	        _normaliser[index] += exp(_exponent[index][yj]);
	      }
	    }
	    for(int yj=0; yj< pow(_Y->rows(),_sizeColValY); yj++)
	    {
	      for(int k=0; k< (*_FM).getMatrixIndexFeat(xi,yj).size();k++)
	      {
	        int index=_FM->getMatrixIndexFeat(xi,yj)[k];
          _expected[index]
					 [(*_FM).getMatrixIndexdX(xi,yj)[k]]
	                 [(*_FM).getMatrixIndexdY(xi,yj)[k]] += exp(_exponent[index][yj])/_normaliser[index];
	      }
	    }
	   }
}

int GISgp:: getIterations()
{
  return _iterations;
}
