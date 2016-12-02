#include "ITMatrix.h"

#include <entropy++/powi.h>

using namespace entropy::iterativescaling;

ITMatrix::ITMatrix(ULContainer *xData,
                   ULContainer *yData,
                   ULContainer *xAlphabet,
                   ULContainer *yAlphabet,
                   ivvector    systX,
                   ivvector    systY,
                   double lambdavalue)
{
  _DataX         = xData;
  _DataY         = yData;
  _xAlphabet    = xAlphabet;
  _yAlphabet    = yAlphabet;
  _systX        = systX;
  _systY        = systY;
  _sizeColDataY = _DataY->columns();
  _sizeColDataX = _DataX->columns();
  _sizeRowDataX = _DataX->rows();
  _sizeRowDataY = _DataY->rows();
  _sizeX        = _xAlphabet->rows();
  _sizeY        = _yAlphabet->rows();
  __featureArray(lambdavalue);

#ifndef MEMORY_EFFICIENT
  __fillX();
  __fillY();
#endif
}

ITMatrix::ITMatrix()
{
  _DataX        = new ULContainer(0,0);
  _DataY        = new ULContainer(0,0);
  _xAlphabet    = new ULContainer(0,0);
  _yAlphabet    = new ULContainer(0,0);
  _sizeColDataY = 0;
  _sizeColDataX = 0;
  _sizeRowDataX = 0;
  _sizeRowDataY = 0;
  _sizeX        = _xAlphabet->rows();
  _sizeY        = _yAlphabet->rows();
  __featureArray(0);
}

ITMatrix::~ITMatrix()
{
  delete _featureArray;
  for(int i=0;i<_systX.size();i++) _systX[i].clear();
  _systX.clear();

  for(int j=0;j<_systY.size();j++) _systY[j].clear();
  _systY.clear();
}

double ITMatrix::getLambda(int i,int ilambdaX, int ilambdaY)
{
  double lambda = _featureArray[i].getLambda(ilambdaX,ilambdaY);
  return lambda;
}

double ITMatrix::getValue(int feat,int rowX,int rowY)
{
  double val=0.0;
  int I = powi(_xAlphabet->rows(),_systX[feat].size());
  int J = powi(_yAlphabet->rows(),_systY[feat].size());
  for(int i = 0; i < I; i++)
  {
    for(int j = 0; j < J; j++)
    {
      val += getDelta(feat,i,j,rowX,rowY)
        * _featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}

double ITMatrix::getValueAlphY(int feat,int rowX,int indexY)
{
  double val = 0.0;
  int I = powi(_xAlphabet->rows(),_systX[feat].size());
  int J = powi(_yAlphabet->rows(),_systY[feat].size());
  for(int i = 0; i < I; i++)
  {
    for(int j = 0; j < J; j++)
    {
      val += getDeltaAlphY(feat,i,j,rowX,indexY)
          * _featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}

double ITMatrix::getValueAlphYAlphX(int feat,int indexX,int indexY)
{
  double val = 0.0;
  int I = powi(_xAlphabet->rows(), _systX[feat].size());
  int J = powi(_yAlphabet->rows(), _systY[feat].size());
  for(int i = 0; i < I; i++)
  {
    for(int j = 0; j < J; j++)
    {
      val += getDeltaAlphYAlphX(feat,i,j,indexX,indexY)
          * _featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}

void ITMatrix::setLambda(int i, int ilambdaX, int ilambdaY,double valuelambda)
{
  _featureArray[i].setLambda(ilambdaX,ilambdaY,valuelambda);
}

int ITMatrix::getDelta(int i,int indexX, int indexY,int rowDataX, int rowDataY)
{
#ifdef MEMORY_EFFICIENT
  int* x = new int[_systX[i].size()];
  int* y = new int[_systY[i].size()];
  index(x, indexX, true,  _systX[i].size());
  index(y, indexY, false,  _systY[i].size());
#endif
  bool equ = true;

  for(int j = 0; j < _systX[i].size(); j++)
  {
#ifdef MEMORY_EFFICIENT
    if((*_DataX)(rowDataX,_systX[i][j]) != x[j])
#else
    if((*_DataX)(rowDataX,_systX[i][j]) != _xFeatureArray[indexX][_systX[i][j]])
#endif
    {
      equ = false;
      break;
    }
  }

  if(equ == true)
  {
    for(int j = 0; j < _systY[i].size(); j++)
    {
#ifdef MEMORY_EFFICIENT
      if((*_DataY)(rowDataY,_systY[i][j]) != y[j]) // _yFeatureArray[indexY][j]
#else
      if((*_DataY)(rowDataY,_systY[i][j]) != _yFeatureArray[indexY][_systY[i][j]])
#endif
      {
        equ = false;
        break;
      }
    }
  }

#ifdef MEMORY_EFFICIENT
  delete[] x;
  delete[] y;
#endif

  if(equ)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

int ITMatrix::getDeltaAlphY(int i, int indexX, int indexY, int rowDataX, int indexDataY)
{
#ifdef MEMORY_EFFICIENT
  int* x = new int[_systX[i].size()];
  int* y = new int[_systY[i].size()];
  int* valy = new int[_sizeColDataY];
  index(x, indexX, true,  _systX[i].size());
  index(y, indexY, false,  _systY[i].size());
  index(valy, indexDataY,false,_sizeColDataY);
#endif

  bool equ = true;
  for(int j = 0; j < _systX[i].size(); j++)
  {
#ifdef MEMORY_EFFICIENT
    if((*_DataX)(rowDataX,_systX[i][j]) != x[j])
#else
    if((*_DataX)(rowDataX,_systX[i][j]) != _xFeatureArray[indexX][_systX[i][j]])
#endif
    {
      equ = false;
      break;
    }
  }

  if(equ == true)
  {
    for(int j=0;j<_systY[i].size();j++)
    {
#ifdef MEMORY_EFFICIENT
      if(valy[_systY[i][j]]!= y[j])
#else
      if(_yFeatureArray[indexDataY][_systY[i][j]] != _yFeatureArray[indexY][_systY[i][j]])
#endif
      {
        equ = false;
        break;
      }
    }
  }

#ifdef MEMORY_EFFICIENT
  delete[] x;
  delete[] y;
  delete[] valy;
#endif
  if(equ)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

int ITMatrix::getDeltaAlphYAlphX(int i, int indexX, int indexY, int indexDataX, int indexDataY)
{
#ifdef MEMORY_EFFICIENT
  int* x    = new int[_systX[i].size()];
  int* y    = new int[_systY[i].size()];
  int* valy = new int[_sizeColDataY];
  int* valx = new int[_sizeColDataX];
  index(x, indexX, true,  _systX[i].size());
  index(y, indexY, false,  _systY[i].size());
  index(valy, indexDataY,false,_sizeColDataY);
  index(valx, indexDataX,true,_sizeColDataX);
#endif

  bool equ = true;

  for(int j = 0; j < _systX[i].size(); j++)
  {
#ifdef MEMORY_EFFICIENT
    if(valx[_systX[i][j]] != x[j])
#else
    if(_xFeatureArray[indexDataX][_systX[i][j]] != _xFeatureArray[indexX][_systX[i][j]])
#endif
    {
      equ = false;
      break;
    }
  }

  if(equ == true)
  {
    for(int j = 0; j < _systY[i].size(); j++)
    {
#ifdef MEMORY_EFFICIENT
      if(valy[_systY[i][j]] != y[j])
#else
      if(_yFeatureArray[indexDataY][_systY[i][j]] != _yFeatureArray[indexY][_systY[i][j]])
#endif
      {
        equ = false;
        break;
      }
    }
  }

#ifdef MEMORY_EFFICIENT
  delete[] x;
  delete[] y;
  delete[] valy;
  delete[] valx;
#endif
  if(equ)
  {
    return 1;
  }
  else
  {
    return -1;
  }
}

// Berechnung der Alphabetreihe aus dem Index
// vector<int> ITMatrix::index(int index, bool x, int sizeCol)
#ifdef MEMORY_EFFICIENT
void ITMatrix::index(int* array, int index, bool x, int sizeCol)
{
  int sizeAlph;
  // vector<int> row;
  double z;
  if(x)
  {
    sizeAlph = _sizeX;
  }
  else
  {
    sizeAlph = _sizeY;
  }

  for(int i = sizeCol; i>0 ; i--)
  {
    z = index % (int)(powi(sizeAlph,i));
    z = z / (powi(sizeAlph,i-1));
    int j = (int)z;
    if(x)
    {
      array[i] = (*_xAlphabet)(j,0);
    }
    else
    {
      array[i] = (*_yAlphabet)(j,0);
    }
  }
  // return row;
}
#endif // MEMORY_EFFICIENT

//ein Array mit den benoetigten Features fuellen
void ITMatrix::__featureArray(double valuelambda)
{
  int sizeX = 0;
  int sizeY = 0;

  _featureArray = new Feature[_systX.size()];

  for(int m=0; m < _systX.size(); m++)
  {
    sizeX = powi(_xAlphabet->rows(), _systX[m].size());
    sizeY = powi(_yAlphabet->rows(), _systY[m].size());
    Feature *K = new Feature(sizeX, sizeY, valuelambda);
    _featureArray[m]= *K;
  }
}


 #ifndef MEMORY_EFFICIENT

void ITMatrix::__fillX()
{
  int nr = powi(_xAlphabet->rows(),_sizeColDataX);
  int nc = _DataX->columns();
  int nx = _sizeColDataX;

  _xFeatureArray = new int*[nr];
  for(int r = 0; r < nr; r++)
  {
    _xFeatureArray[r] = new int[nc];
  }

  for(int r = 0; r < nr; r++)
  {
    for(int c = 0; c < nc; c++)
    {
      float f  = powi(nx, c);
      int   rf = (int)(((float)r) / f);
      int   xindex = rf % nx;
      _xFeatureArray[r][nc-c-1] = _xAlphabet->get(xindex, 0);
    }
  }
}

void ITMatrix::__fillY()
{
  int nr = powi(_yAlphabet->rows(), _sizeColDataY);
  int nc = _sizeColDataY;
  int nx = _yAlphabet->rows();

  _yFeatureArray = new int*[nr];
  for(int r = 0; r < nr; r++)
  {
    _yFeatureArray[r] = new int[nc];
  }

  for(int r = 0; r < nr; r++)
  {
    for(int c = 0; c < nc; c++)
    {
      float f  = powi(nx, c);
      int   rf = (int)(((float)r) / f);
      int   xindex = rf % nx;
      _yFeatureArray[r][nc-c-1] = _yAlphabet->get(xindex, 0);
    }
  }
}

 #endif // MEMORY_EFFICIENT
