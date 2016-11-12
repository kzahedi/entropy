#include "GISgp.h"

using namespace entropy::iterativescaling;

//training data, alphabete , startwert fuer lambda, startwert fuer delta, wert fuer sigma, test auf time, sekunden fuer den test
GISgp::GISgp(DContainer &xData,
             DContainer &yData,
             DContainer &xAlphabet,
             DContainer &yAlphabet,
             ivvector systX,
             ivvector systY,
             IsParameter param)
  : IterativeScalingBase(xData, yData, xAlphabet, yAlphabet, systX, systY, param, true)
{
  _fm         = (FeatureMatrix*)_imatrix;
  _normaliser = new double[_sizeSystX];
  _exponent   = new double*[_sizeSystX];

  int Y = (int)pow(_yAlphabet->rows(),_sizeColDataY);
  for(int i = 0; i < _sizeSystX; i++)
  {
    _exponent[i]= new double[Y];
  }

  _expected = new double**[_sizeSystX];
  for(int i = 0; i < _sizeSystX; i++)
  {
    int K = (int) pow(_xAlphabet->rows(),_systX[i].size());
    int L = (int) pow(_yAlphabet->rows(),_systY[i].size());
    _expected[i] = new double*[K];
    for(int k = 0; k < K; k++)
    {
      _expected[i][k] = new double[L];
      for(int l = 0; l < _sizeY; l++)
      {
        _expected[i][k][l] = 0;
      }
    }
  }

  _delta = new double**[_sizeSystX];
  for(int i = 0; i < _sizeSystX; i++)
  {
    int J = (int) pow(_xAlphabet->rows(),_systX[i].size());
    int K = (int) pow(_yAlphabet->rows(),_systY[i].size());
    _delta[i] = new double*[J];
    for(int j = 0; j < J; j++)
    {
      _delta[i][j]= new double[K];
      for(int k = 0; k < K; k++)
      {
        _delta[i][j][k] = param.lambdadeltaval;
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

GISgp::GISgp(ULContainer &xData,
             ULContainer &yData,
             DContainer &xAlphabet,
             DContainer &yAlphabet,
             ivvector systX,
             ivvector systY,
             IsParameter param)
  : IterativeScalingBase(xData, yData, xAlphabet, yAlphabet, systX, systY, param, true)
{
  _fm         = (FeatureMatrix*)_imatrix;
  _normaliser = new double[_sizeSystX];
  _exponent   = new double*[_sizeSystX];

  int N = (int) pow(_yAlphabet->rows(),_sizeColDataY);
  for(int i = 0; i < _sizeSystX; i++)
  {
    _exponent[i] = new double[N];
  }

  _expected = new double**[_sizeSystX];
  for(int i = 0; i < _sizeSystX; i++)
  {
    int K = (int)pow(_xAlphabet->rows(),_systX[i].size());
    _expected[i] = new double*[K];
    for(int k = 0; k < K; k++)
    {
      int L = (int)pow(_yAlphabet->rows(),_systY[i].size());
      _expected[i][k] = new double[L];
      for(int l = 0; l < _sizeY; l++)
      {
        _expected[i][k][l] = 0;
      }
    }
  }

  _delta = new double**[_sizeSystX];
  for(int i = 0; i < _sizeSystX; i++)
  {
    int J = (int) pow(_xAlphabet->rows(),_systX[i].size());
    int K = (int) pow(_yAlphabet->rows(),_systY[i].size());
    _delta[i] = new double*[J];
    for(int j = 0; j < J; j++)
    {
      _delta[i][j] = new double[K];
      for(int k = 0; k < K; k++)
      {
        _delta[i][j][k] = param.lambdadeltaval;
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
    for(int i = 0; i < _sizeSystX; i++)
    {
      int K = pow(_xAlphabet->rows(),_systX[i].size());
      for(int k = 0; k < K; k++)
      {
        delete[] _expected[i][k];
      }
      delete[] _expected[i];
    }
    delete[] _expected;
  }

  if(_delta != NULL)
  {
    for(int i = 0; i < _sizeSystX; i++)
    {
      for(int j = 0; j < pow(_xAlphabet->rows(),_systX[i].size()); j++)
      {
        delete[] _delta[i][j];
      }
      delete[] _delta[i];
    }
    delete[] _delta;
  }

  delete [] _exponent;
  delete [] _normaliser;

  _conv.clear();
}

double GISgp::__calculateIteration(double featconst, double sigma, bool test)
{
  double l = 0;
  __getexp();
  for(int feat = 0; feat < _sizeSystX; feat++)
  {
    int DI = pow(_xAlphabet->rows(),_systX[feat].size());
    int DJ = pow(_yAlphabet->rows(),_systY[feat].size());
    for(int deltai = 0; deltai < DI; deltai++)
    {
      for(int deltaj = 0; deltaj < DJ; deltaj++)
      {
        double newl = 0.0;
        double oldl = _fm->getFeatureArraylambda(feat,deltai,deltaj);
        double zOld = 2;
        double z    = 1;
        while(fabs(z-zOld) > 0.0001)
        {
          zOld     = z;
          z        = (oldl+_delta[feat][deltai][deltaj])/pow(sigma,2)
            + _expected[feat][deltai][deltaj]*exp(_delta[feat][deltai][deltaj]*featconst)
            - _observed[feat][deltai][deltaj];
          double n = + 1/(pow(sigma,2)) // TODO = + ?
            + _expected[feat][deltai][deltaj]
            * featconst
            * exp(_delta[feat][deltai][deltaj]* featconst);

            if(fabs(n)<0.00000001){ cout << "  hieer " << endl;} // TODO Exeption?
            _delta[feat][deltai][deltaj] = _delta[feat][deltai][deltaj]-(z/n);
        }

        newl = oldl + _delta[feat][deltai][deltaj];
        _fm->setFeatureArraylambda(feat,deltai,deltaj,newl);
        l += fabs((oldl+_delta[feat][deltai][deltaj])
                  / pow(sigma,2)
                  + _expected[feat][deltai][deltaj]
                  * exp(_delta[feat][deltai][deltaj]*featconst)
                  - _observed[feat][deltai][deltaj]);
      }
    }
  }
  _iterations++;
  if(test) _conv.push_back(l);
  return l;
}

void GISgp::__gisgp(int maxit, double konv, double sigma, bool test,int seconds)
{
  //constant c for delta
  double featconst = __getFeatconst();

  int N = pow(_yAlphabet->rows(),_sizeColDataY);
  for(int i = 0; i < _sizeSystX; i++)
  {
    for(int j = 0; j < N; j++)
    {
      _exponent[i][j]= 0.0;
    }
  }
  double l     = 1;
  double utime = 0;
  time_t befor;
  time_t after;
  _iterations = 0;
  while(utime < seconds)
  {
    befor  = time(NULL);
    __calculateIteration(featconst,sigma,test);
    after  = time(NULL);
    utime += difftime(after,befor);
  }
}

void GISgp::__gisgp(int maxit, double konv, double sigma, bool test)
{
  //constant c for delta
  double featconst = __getFeatconst();
  int    N         = pow(_yAlphabet->rows(),_sizeColDataY);

  for(int i = 0; i < _sizeSystX; i++)
  {
    for(int j = 0; j < N; j++)
    {
      _exponent[i][j] = 0.0;
    }
  }

  double l = konv + 1.0;
  _iterations = 0;
  while(_iterations < maxit && fabs(l) >= konv)
  {
    l = __calculateIteration(featconst,sigma,true);
  }
}

double GISgp::getconv(int i)
{
  return _conv[i];
}

int GISgp::getsizeconv()
{
  return _conv.size();
}

double GISgp::__getFeatconst()
{
  double r = 0.0;
  int N = pow(_yAlphabet->rows(),_sizeColDataY);
  for(int i = 0; i < _sizeRowDataX; i++) // i-th data row
  {
    for(int j = 0; j < N; j++) // y-alphabet
  {
    int v = _fm->getMatrixIndexFeat(i,j).size(); // the number of matching deltas
    if(v > r) r = v;
  }
  }
  return r;
}

void GISgp::__getexp()
{
  for(int i=0; i<_sizeSystX; i++)
  {
    int K = (int) pow(_xAlphabet->rows(),_systX[i].size());
    int L = (int) pow(_yAlphabet->rows(),_systY[i].size());
    for(int k = 0; k < K; k++)
    {
      for(int l = 0; l < L;l++)
      {
        _expected[i][k][l] = 0;
      }
    }
  }

  int YJ = pow(_yAlphabet->rows(),_sizeColDataY);
  for(int xi = 0; xi < _sizeRowDataX; xi++)
  {
    for(int i = 0; i < _sizeSystX; i++)
    {
      _normaliser[i] = 0;
    }
    for(int yj = 0; yj < YJ; yj++)
    {
      for(int k = 0; k < _fm->getMatrixIndexFeat(xi,yj).size(); k++)
      {
        int index = _fm->getMatrixIndexFeat(xi,yj)[k];
        _exponent[index][yj] =_fm->getFeatureArrayvalueAlphY(index,xi,yj);
        _normaliser[index]  += exp(_exponent[index][yj]);
      }
    }

    for(int yj = 0; yj < YJ; yj++)
    {
      for(int k=0; k< _fm->getMatrixIndexFeat(xi,yj).size();k++)
      {
        int index=_fm->getMatrixIndexFeat(xi,yj)[k];
        _expected[index]
          [_fm->getMatrixIndexdX(xi,yj)[k]]
          [_fm->getMatrixIndexdY(xi,yj)[k]] += exp(_exponent[index][yj])/_normaliser[index];
      }
    }
  }
}

int GISgp::getIterations()
{
  return _iterations;
}
