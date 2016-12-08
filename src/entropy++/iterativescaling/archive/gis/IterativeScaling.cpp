#include "IterativeScaling.h"

#include <entropy++/powi.h>

using namespace entropy::iterativescaling::gis;

#define EPSILON 0.00000001


IterativeScaling::IterativeScaling(ULContainer *xData,
         ULContainer *yData,
         ULContainer *xAlphabet,
         ULContainer *yAlphabet,
         ivvector systX,
         ivvector systY,
         IsParameter param)
  : IterativeScalingBase(xData, yData, xAlphabet, yAlphabet, systX, systY, param, true)
{
  _param    = param;
  _expected = new double**[_sizeSystX];
  _fm       = (FeatureMatrix*)_imatrix;
  int K = 0;
  int L = 0;
  for(int i=0; i<_sizeSystX; i++)
  {
    K = (int)powi(_xAlphabet->rows(), _systX[i].size());
    L = (int)powi(_yAlphabet->rows(), _systY[i].size());
    _expected[i] = new double*[K];
    for(int k=0; k < K; k++)
    {
      _expected[i][k]= new double[L];
      for(int l = 0; l < L; l++)
      {
        _expected[i][k][l] = 0;
      }
    }
  }
  _normaliser = new double[_sizeSystX];
  _exponent   = new double*[_sizeSystX];
  int size = (int)powi(_yAlphabet->rows(), _sizeColDataY);
  for(int i = 0; i < _sizeSystX; i++)
  {
    _exponent[i] = new double[size];
  }

  // cout << "Data X:" << endl << xData << endl << "Data Y: " << endl << yData << endl;
  if(param.konvtime)
  {
    __gis( param.konv, param.seconds, param.test);
  }
  else
  {
    if(param.time ) __gis(param.seconds, param.test);
    else            __gis(param.maxit, param.konv, param.test);
  }
}

IterativeScaling::~IterativeScaling()
{
  if(_expected != NULL)
  {
    for(int i = 0; i < _sizeSystX; i++)
    {
      int K = powi(_xAlphabet->rows(), _systX[i].size());
      for(int k = 0; k < K; k++)
      {
        delete[] _expected[i][k];
      }
      delete[] _expected[i];
    }
    delete[] _expected;
  }

  if(_exponent != NULL)
  {
    for(int i = 0; i < _sizeSystX; i++){
      delete[] _exponent[i];
    }
    delete[] _exponent;
  }
  delete[] _normaliser;
}
//konstante zur Verringerung der Schrittgroesse
double IterativeScaling::__getFeatconst()
{
  double r = 0.0;
  int N = powi(_yAlphabet->rows(), _sizeColDataY);
  int sizeUniqueX= _fm->getSizeUnique();
  for(int i = 0; i < sizeUniqueX; i++) // i-th data row
  {
    for(int j = 0; j < N; j++) // y-alphabet
    {
      int v = _fm->getMatrixIndexFeatSize(i,j); // the number of matching deltas
      if(v > r) r = v;
    }
  }
  return r;
}
// Mit den derzeitigen Parametern lambda berechnete Anzahl der Vorkommen von feature(deltai,deltaj)
void IterativeScaling::__getExpected()
{
  for(int i = 0; i < _sizeSystX; i++)
  {
    int K = (int)powi(_xAlphabet->rows(), _systX[i].size());
    int L = (int)powi(_yAlphabet->rows(), _systY[i].size());
    for(int k = 0; k < K; k++)
    {
      for(int l = 0; l < L; l++)
      {
        _expected[i][k][l] = 0;
      }
    }
  }
  int featsize;
  int indexUniqueX;
  int YJ = powi(_yAlphabet->rows(), _sizeColDataY);
  for(int xi = 0; xi < _sizeRowDataX; xi++)
  {
    indexUniqueX= _fm->getUniqueIndex(xi); // index der Reihe aus UniqueX, die mit xi aus DataX uebereinstimmt
    for(int i = 0; i < _sizeSystX; i++)
    {
      _normaliser[i] = 0.0;
    }
    for(int yj = 0; yj < YJ; yj++)
    {
      featsize = _fm->getMatrixIndexFeatSize(indexUniqueX,yj);
      for(int k = 0; k < featsize; k++)
      {
        int index = _fm->getMatrixIndexFeatValue(indexUniqueX,yj,k);
        _exponent[index][yj] = _fm->getValueAlphY(index,xi,yj);
        _normaliser[index]  +=  exp(_exponent[index][yj]);
      }
    }
    for(int yj = 0; yj < YJ; yj++)
    {
      featsize =  _fm->getMatrixIndexFeatSize(indexUniqueX, yj);
      for(int k = 0; k < featsize ;k++)
      {
        int index = _fm->getMatrixIndexFeatValue(indexUniqueX, yj, k);
        _expected[index]
          [_fm->getMatrixIndexdXValue(indexUniqueX,yj,k)]
          [_fm->getMatrixIndexdYValue(indexUniqueX,yj,k)] += exp(_exponent[index][yj])/_normaliser[index];
      }
    }
  }
}
// auf zeit
void IterativeScaling::__gis(int seconds, bool test) // TODO remove test
{
  //constant c for delta
  double featconst = __getFeatconst();

  for(int k = 0; k < _sizeSystX; k++)
  {
    int I = powi(_yAlphabet->rows(), _sizeColDataY);
    for(int i = 0; i < I; i++)
    {
      _exponent[k][i] = 0.0;
    }
  }
  _normaliser  = new double[_sizeSystX];
  _iterations  = 0;
  double l     = 1;
  double utime = 0;
  time_t befor;
  time_t after;
  while(utime < seconds)
  {
    befor  = time(NULL);
    __calculateIteration(featconst,test);
    after  = time(NULL);
    utime += difftime(after,befor);
  }
}

// auf zeit und konv
void IterativeScaling::__gis( double konv,int seconds, bool test)
{
  double featconst = __getFeatconst();

  int size = powi(_yAlphabet->rows(), _sizeColDataY);
  for(int k = 0; k < _sizeSystX; k++)
  {
    for(int i = 0; i < size; i++)
    {
      _exponent[k][i]=0.0;
    }
  }
  _normaliser  = new double[_sizeSystX];
  _iterations  = 0;
  double l     = konv + 1.0;
  double utime = 0;
  time_t befor;
  time_t after;
  while(utime < seconds && fabs(l) >= konv)
  {
    befor  = time(NULL);
    l      = __calculateIteration(featconst,test);
    after  = time(NULL);
    utime += difftime(after,befor);
  }
}

// nach iterationen und konv
void IterativeScaling::__gis(int maxit, double konv,bool test)
{
  //constant c for delta
  double featconst = __getFeatconst();

  for(int k = 0; k < _systX.size(); k++)
  {
    int I = powi(_yAlphabet->rows(),_sizeColDataY);
    _exponent[k] = new double[I]; // lambda_i * f_i
    for(int i = 0; i < I; i++)
    {
      _exponent[k][i] = 0.0;
    }
  }
  double l    = konv + 1.0;
  _iterations = 0;
  while(_iterations < maxit && fabs(l) >= konv)
  {
    cout << "iteration " << _iterations << endl;
    l = __calculateIteration(featconst, test);
    cout << l << endl;
  }
} 

double IterativeScaling::__calculateIteration(double featconst, bool test)
{
  double l = 0;
  double oldl = 0.0;
  double newl = 0.0;
  __getExpected();
  for(int feat = 0; feat < _sizeSystX; feat++)
  {
    // jedes delta hat ein x_i und ein y_j
    int DI = powi(_xAlphabet->rows(),_systX[feat].size());
    int DJ = powi(_yAlphabet->rows(),_systY[feat].size());
    for(int deltai = 0; deltai < DI; deltai++)
    {
      for(int deltaj = 0; deltaj < DJ; deltaj++)
      {
        newl = 0.0;
        if(fabs(_expected[feat][deltai][deltaj]) < EPSILON &&
           fabs(_observed[feat][deltai][deltaj]) > EPSILON)
        {
          _expected[feat][deltai][deltaj] = 0.01; // TODO check if other values might be better
        }
        if(fabs(_observed[feat][deltai][deltaj]) > EPSILON)
        {
          // TODO 0.1 as learning rate parameter
          oldl = _fm->getLambda(feat, deltai, deltaj);
          newl = oldl
            + 0.1*(1.0/featconst)
            * log(_observed[feat][deltai][deltaj]/_expected[feat][deltai][deltaj]);
          _fm->setLambda(feat,deltai,deltaj,newl);
        }
        l += fabs(_observed[feat][deltai][deltaj] - _expected[feat][deltai][deltaj]);
      }
    }
  }

  _iterations++;
  if(test) _conv.push_back(l);
  return l;
}

int IterativeScaling::getIterations()
{
  return _iterations;
}

double IterativeScaling::getconv(int i)
{
  return _conv[i];
}

int IterativeScaling::getsizeconv()
{
  return _conv.size();
}

