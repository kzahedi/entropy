#include "SCGIS.h"

using namespace entropy::iterativescaling;

#define EPSILON 0.00000001

SCGIS::SCGIS(DContainer &xData,
             DContainer &yData,
             DContainer &xAlphabet,
             DContainer &yAlphabet,
             ivvector systX,
             ivvector systY,
             IsParameter param)
 : IterativeScalingBase(xData, yData, xAlphabet, yAlphabet, systX, systY, param, false)
{
  _param    = param;
  _exponent = new double**[_sizeSystX];
  _im       = (InstanceMatrix*)_imatrix;
  int Y = (int)pow(_yAlphabet->rows(),_sizeColDataY);

  for(int i = 0; i < _sizeSystX; i++)
  {
    _exponent[i] = new double*[_sizeRowDataX];
    for(int xi = 0; xi < _sizeRowDataX; xi++)
    {
      _exponent[i][xi] = new double[Y];
      for(int y = 0; y < Y; y++)
      {
        _exponent[i][xi][y] = 0;
      }
    }
  }

  _normaliser = new double*[_sizeSystX];
  for(int i = 0; i < _sizeSystX; i++)
  {
    _normaliser[i] = new double[_sizeRowDataX];
    for(int k = 0; k < _sizeRowDataX;k++)
    {
      _normaliser[i][k] = Y;
    }
  }

  _delta=0.0;
  if(param.konvtime)
  {
     __scgis(param.konv, param.seconds, param.test);
  }
  else{
    if(param.time )
    {
       __scgis(param.seconds, param.test);
    }
    else
    {
     __scgis(param.maxit, param.konv, param.test);
    }
  }
}

SCGIS::SCGIS(ULContainer &xData,
             ULContainer &yData,
             DContainer &xAlphabet,
             DContainer &yAlphabet,
             ivvector systX,
             ivvector systY,
             IsParameter param)
  : IterativeScalingBase(xData, yData, xAlphabet, yAlphabet, systX, systY, param, false)
{
  _im       = (InstanceMatrix*)_imatrix;
  _param    = param;
  _exponent = new double**[_sizeSystX];
  int Y     = (int)pow(_yAlphabet->rows(),_sizeColDataY);
  for(int i = 0; i < _sizeSystX; i++)
  {
    _exponent[i] = new double*[_sizeRowDataX];
    for(int xi = 0; xi < _sizeRowDataX; xi++)
    {
      _exponent[i][xi] = new double[Y];
      for(int y = 0; y < Y; y++)
      {
        _exponent[i][xi][y]=0;
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

  _delta=0.0;
  if(param.konvtime)
  {
    __scgis(param.konv, param.seconds, param.test); // TODO rename functions
  }
  else{
    if(param.time) __scgis(param.seconds, param.test);
    else           __scgis(param.maxit, param.konv, param.test);
  }
}


double SCGIS::getconv(int i)
{
  return _conv[i];
}

int SCGIS::getsizeconv()
{
  return _conv.size();
}

SCGIS::~SCGIS()
{
  for(int i=0;i<_sizeSystX;i++)
  {
    for(int j=0;j<_sizeRowDataX;j++)
    {
      delete[] _exponent[i][j];
    }
    delete[] _exponent[i];
  }
  delete[] _exponent;

  for(int i=0;i<_sizeSystX;i++){
    delete[] _normaliser[i];
  }
  delete[] _normaliser;

  _conv.clear();
}


void SCGIS::__scgis( double konv,int seconds, bool test)
{
  double l     = 1;
  double utime = 0;
  _iterations  = 0;
  time_t befor;
  time_t after;

  while(utime < seconds && fabs(l) > konv)
  {
    befor  = time(NULL);
    l      = __calculateIteration(test);
    after  = time(NULL);
    utime += difftime(after,befor);
  }
}

void SCGIS::__scgis(int seconds, bool test)
{
  double l     = 1;
  double utime = 0;
  _iterations  = 0;
  time_t befor;
  time_t after;

  while(utime<seconds)
  {
    befor  = time(NULL);
    __calculateIteration(test);
    after  = time(NULL);
    utime += difftime(after,befor);
  }

}
void SCGIS::__scgis(int maxit, double konv, bool test)
{
  double l    = konv + 1.0;
  _iterations = 0;
  while(_iterations < maxit && fabs(l) > konv) 
  {
    l = __calculateIteration(test);
  }
}

double SCGIS::__calculateIteration(bool test)
{
  double l = 0.0;
  int    Y = pow(_sizeY,_sizeColDataY);
  for(int feat = 0; feat < _sizeSystX; feat++)
  {
    int DI = pow(_sizeX,_systX[feat].size());
    int DJ = pow(_sizeY,_systY[feat].size());
    for(int delti = 0; delti < DI; delti++)
    {
      for(int deltj = 0; deltj < DJ; deltj++)
      {
        double expected = 0.0;
        for(int y = 0; y < Y; y++)
        {
          for(int k = 0; k < _im->getInstanceMatrixX(feat,delti,deltj).size(); k++)
          {
            if(_im->getInstanceMatrixY(feat,delti,deltj)[k] == y)
            {
              int x = _im->getInstanceMatrixX(feat,delti,deltj)[k];
              if(fabs(_normaliser[feat][x]) > EPSILON )
              {
                // f_i(\bar{x}_j, y) is given by the 2 for and 1 if above
                expected += exp(_exponent[feat][x][y])/_normaliser[feat][x];
              }
            }
          }
        }

        double newl = 0.0;
        if(fabs(expected) < EPSILON &&
           _observed[feat][delti][deltj] > EPSILON)
        {
          expected=0.01;
        }

        if(fabs(_observed[feat][delti][deltj]) > EPSILON)
        {
          _delta = log(_observed[feat][delti][deltj]/expected);
          newl   = _im->getFeatureArraylambda(feat,delti,deltj) + _delta;
          _im->setFeatureArraylambda(feat,delti,deltj,newl);
        }
        else
        {
          _delta = -1; // TODO evtl besseren Wert finden
        }

        l += fabs(_observed[feat][delti][deltj] - expected);
        for(int y = 0; y < Y; y++)
        {
          for(int k = 0; k < _im->getInstanceMatrixX(feat,delti,deltj).size(); k++)
          {
            if(_im->getInstanceMatrixY(feat,delti,deltj)[k] == y)
            {
              int x = _im->getInstanceMatrixX(feat,delti,deltj)[k];
              _normaliser[feat][x]  -= exp(_exponent[feat][x][y]);
              _exponent[feat][x][y] += _delta;
              _normaliser[feat][x]  += exp(_exponent[feat][x][y]);
            }
          }
        }
      }
    }
  }

  _iterations++;
  if(test) _conv.push_back(l);

  //  cout << l << " " << _iterations << endl;
  return l;
}

int SCGIS::getIterations()
{
  return _iterations;
}

