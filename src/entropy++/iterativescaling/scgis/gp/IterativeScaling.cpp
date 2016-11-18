#include "IterativeScaling.h"

#include <entropy++/powi.h>

using namespace entropy::iterativescaling::scgis::gp;


#define EPSILON 0.00000001

IterativeScaling::IterativeScaling(ULContainer *xData,
                 ULContainer *yData,
                 ULContainer *xAlphabet,
                 ULContainer *yAlphabet,
                 ivvector systX,
                 ivvector systY,
                 IsParameter param)
  : IterativeScalingBase(xData, yData, xAlphabet, yAlphabet, systX, systY, param, false)
{
  _im       = (InstanceMatrix*)_imatrix;
  _exponent = new double**[_sizeSystX];
  int Y     = powi(_yAlphabet->rows(),_sizeColDataY);
  for(int i = 0; i < _sizeSystX; i++)
  {
    _exponent[i] = new double*[_sizeRowDataX];
    for(int xi = 0; xi < _sizeRowDataX; xi++)
    {
      _exponent[i][xi] = new double[Y];
      for(int y = 0; y < Y; y++)
      {
        _exponent[i][xi][y] = 0.0;
      }
    }
  }
  _normaliser = new double*[_sizeSystX];
  for(int i = 0; i < _sizeSystX; i++)
  {
    _normaliser[i] = new double[_sizeRowDataX];
    for(int k = 0; k < _sizeRowDataX; k++)
    {
      _normaliser[i][k] = Y;
    }
  }
  _delta = new double**[_sizeSystX];
  for(int i = 0; i < _sizeSystX; i++)
  {
    int K = (int)powi(_xAlphabet->rows(),_systX[i].size());
    int L = (int)powi(_yAlphabet->rows(),_systY[i].size());
    _delta[i] = new double*[K];
    for(int k = 0; k < K; k++)
    {
      _delta[i][k] = new double[L];
      for(int l = 0; l < L; l++)
      {
        _delta[i][k][l] = param.lambdadeltaval;
      }
    }
  }
  // TODO rename functions
  if(param.time) __scgis(param.maxit,param.konv,param.test,param.sigma,param.seconds);
  else           __scgis(param.maxit,param.konv,param.test,param.sigma);
}

IterativeScaling::~IterativeScaling()
{
  int J = powi(_xAlphabet->rows(),_sizeColDataX);
  for(int i = 0; i < _sizeSystX; i++)
  {
    for(int j = 0; j < J; j++)
    {
      delete[] _exponent[i][j];
    }
    delete[] _exponent[i];
  }
  delete[] _exponent;

  for(int i = 0; i < _sizeSystX; i++)
  {
    delete[] _normaliser[i];
  }
  delete[] _normaliser;

  for(int i = 0; i < _sizeSystX; i++)
  {
    int J = powi(_xAlphabet->rows(),_systX[i].size());
    for(int j = 0; j < J; j++)
    {
      delete[] _delta[i][j];
    }
    delete[] _delta[i];
  }
  delete[] _delta;
  _conv.clear();
}

double IterativeScaling::getconv(int i)
{
  return _conv[i];
}

int IterativeScaling::getsizeconv()
{
  return _conv.size();
}

double IterativeScaling::__calculateIteration(bool test, double sigma)
{
  double l = 0.0;
  int Y = powi(_yAlphabet->rows(),_sizeColDataY);
  for(int feat = 0; feat < _sizeSystX; feat++)
  {
    int DI = powi(_xAlphabet->rows(),_systX[feat].size());
    int DJ = powi(_yAlphabet->rows(), _systY[feat].size());
    for(int delti = 0; delti < DI; delti++)
    {
      for(int deltj = 0; deltj < DJ; deltj++)
      {
        double expected=0;
        for(int y = 0; y < Y; y++)
        {
          for(int k = 0; k < _im->getInstanceMatrixX(feat,delti,deltj).size(); k++)
          {
            if(_im->getInstanceMatrixY(feat,delti,deltj)[k] == y)
            {
              int x = _im->getInstanceMatrixX(feat,delti,deltj)[k];
              if(fabs(_normaliser[feat][x]) > EPSILON)
              {
                expected += exp(_exponent[feat][x][y])/_normaliser[feat][x];
              }
            }
          }
        }
        double newl = 0.0;
        double oldl = _im->getFeatureArraylambda(feat,delti,deltj);
        double zOld = 2;
        double z    = 1;
        while(fabs(z-zOld) > 0.0001)
        {
          zOld = z;
          z    = (oldl+_delta[feat][delti][deltj])/powi(sigma,2)
            + expected*exp(_delta[feat][delti][deltj]) - _observed[feat][delti][deltj];
          double n = 1.0/(powi(sigma,2)) + expected*exp(_delta[feat][delti][deltj]);
          _delta[feat][delti][deltj] = _delta[feat][delti][deltj] - (z/n);
        }
        newl  = _im->getFeatureArraylambda(feat,delti,deltj)+_delta[feat][delti][deltj];
        l    += fabs((oldl+_delta[feat][delti][deltj])/powi(sigma,2)
                     + expected*exp(_delta[feat][delti][deltj])
                     - _observed[feat][delti][deltj]);
        _im->setFeatureArraylambda(feat,delti,deltj,newl);
        for(int y = 0; y < powi(_yAlphabet->rows(),_sizeColDataY); y++)
        {
          for(int k = 0; k < _im->getInstanceMatrixX(feat,delti,deltj).size(); k++)
          {
            if(_im->getInstanceMatrixY(feat,delti,deltj)[k] == y)
            {
              int x = _im->getInstanceMatrixX(feat,delti,deltj)[k];
              _normaliser[feat][x]  -= exp(_exponent[feat][x][y]);
              _exponent[feat][x][y] += _delta[feat][delti][deltj];
              _normaliser[feat][x]  += exp(_exponent[feat][x][y]);
            }
          }
        }
      }
    }
  }
  _iterations++;
  if(test) _conv.push_back(l);
  return l;
}

void IterativeScaling::__scgis(int maxit, double konv, bool test, double sigma, int seconds)
{
  double utime = 0;
  _iterations  = 0;
  time_t befor;
  time_t after;
  while(utime < seconds)
  {
    befor  = time(NULL);
    __calculateIteration(test,sigma);
    after  = time(NULL);
    utime += difftime(after,befor);
  }
}

void IterativeScaling::__scgis(int maxit, double konv, bool test, double sigma)
{
  double l    = 1;
  _iterations = 0;
  while(_iterations < maxit) //&& fabs(l)>konv
  {
    l = __calculateIteration(test,sigma);
    // cout << l << endl; // TODO required?
  }
}

int IterativeScaling::getIterations()
{
  return _iterations;
}
