#include "ITMatrix.h"

#include <entropy++/powi.h>

using namespace entropy::iterativescaling;


ITMatrix::ITMatrix(ULContainer *xData,
                   ULContainer *yData,
                   ULContainer *xAlphabet,
                   ULContainer *yAlphabet,
                   vector<vector<int> > systX,
                   vector<vector<int> > systY,
                   double lambdavalue)
{
  // assert(systX.size() == systY.size());
  _valX         = xData;
  _valY         = yData;
  _xAlphabet    = xAlphabet;
  _yAlphabet    = yAlphabet;
  _cmi          = true;
  _systX        = systX;
  _systY        = systY;
  _sizeColDataY = _valY->columns();
  _sizeColDataX = _valX->columns();
  _sizeRowDataX = _valX->rows();
  _sizeRowDataY = _valY->rows();
  _sizeX        = _xAlphabet->rows();
  _sizeY        = _yAlphabet->rows();
  __featureArray(lambdavalue);
}

ITMatrix::ITMatrix()
{
  _valX         = new ULContainer(0,0);
  _valY         = new ULContainer(0,0);
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

double ITMatrix::getFeatureArraylambda(int i,int ilambdaX, int ilambdaY)
{
  // assert(i<_systX.size());  //ilambdaX und ilambdaY werden in Feature geprueft
  double lambda = _featureArray[i].getLambda(ilambdaX,ilambdaY);
  return lambda;
}

double ITMatrix::getFeatureArrayvalue(int feat,int rowX,int rowY)
{
  // assert(feat<_systX.size());
  double val=0.0;
  int I = powi(_xAlphabet->rows(),_systX[feat].size());
  int J = powi(_yAlphabet->rows(),_systY[feat].size());
  for(int i = 0; i < I; i++)
  {
    for(int j = 0; j < J; j++)
    {
      val += getFeatureArraydelta(feat,i,j,rowX,rowY)
        * _featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}

double ITMatrix::getFeatureArrayvalueAlphY(int feat,int rowX,int indexY)
{
  // assert(feat<_systX.size());
  double val = 0.0;
  int I = powi(_xAlphabet->rows(),_systX[feat].size());
  int J = powi(_yAlphabet->rows(),_systY[feat].size());
  for(int i = 0; i < I; i++)
  {
    for(int j = 0; j < J; j++)
    {
      val += getFeatureArraydeltaAlphY(feat,i,j,rowX,indexY)
          * _featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}

double ITMatrix::getFeatureArrayvalueAlphYAlphX(int feat,int indexX,int indexY)
{
  // assert(feat<_systX.size());
  double val = 0.0;
  int I = powi(_xAlphabet->rows(), _systX[feat].size());
  int J = powi(_yAlphabet->rows(), _systY[feat].size());
  for(int i = 0; i < I; i++)
  {
    for(int j = 0; j < J; j++)
    {
      val += getFeatureArraydeltaAlphYAlphX(feat,i,j,indexX,indexY)
          * _featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}

void ITMatrix::setFeatureArraylambda(int i, int ilambdaX, int ilambdaY,double valuelambda)
{
  // assert(i<_systX.size());
  _featureArray[i].setLambda(ilambdaX,ilambdaY,valuelambda);
}

int ITMatrix::getFeatureArraydelta(int i,int indexX, int indexY,int rowDataX, int rowDataY)
{
  // assert(i<_systX.size());
  // assert(indexX < pow(_xAlphabet->rows(),_systX[i].size()));
  // assert(indexY < pow(_yAlphabet->rows(),_systY[i].size()));
  // assert(rowDataX < _sizeRowDataX);
  // assert(rowDataY < _sizeRowDataY);
  double* x = new double[_systX[i].size()];
  double* y = new double[_systY[i].size()];
  index(x, indexX, true,  _systX[i].size());
  index(y, indexY, false,  _systY[i].size());
  // vector<int> x = index(indexX, true,  _systX[i].size());
  // vector<int> y = index(indexY, false, _systY[i].size());
  bool equ = true;

  for(int j = 0; j < _systX[i].size(); j++)
  {
    // assert(_systX[i][j] <= _sizeColDataX);
    if(_cmi == false)
    {
      if((*_valX)(rowDataX,_systX[i][j]) != x[j])
      {
        equ = false;
      }
    }
    else
    {
      if((*_valX)(rowDataX,_systX[i][j]) != x[j])
      {
        equ = false;
      }
    }
  }

  for(int j = 0; j < _systY[i].size(); j++)
  {
    // assert(_systY[i][j] <= _sizeColDataY);
    if(_cmi == false)
    {
      if((*_valY)(rowDataY,_systY[i][j]) != y[j])
      {
        equ = false;
      }
    }
    else
    {
      if((*_valY)(rowDataY,_systY[i][j]) != y[j])
      {
        equ = false;
      }
    }
  }

  if(equ)
  {
    return 1;
  }
  else
  {
    return -1;
  }
  delete[] x;
  delete[] y;
}

int ITMatrix::getFeatureArraydeltaAlphY(int i, int indexX, int indexY, int rowDataX, int indexDataY)
{
  // assert(i<_systX.size());
  // assert(indexX < pow(_xAlphabet->rows(),_systX[i].size()));
  // assert(indexY < pow(_yAlphabet->rows(),_systY[i].size()));
  // assert(rowDataX < _sizeRowDataX);
  // assert(indexDataY < pow(_yAlphabet->rows(),_sizeColDataY));
  double* x = new double[_systX[i].size()];
  double* y = new double[_systY[i].size()];
  double* valy = new double[_sizeColDataY];
  index(x, indexX, true,  _systX[i].size());
  index(y, indexY, false,  _systY[i].size());
  index(valy, indexDataY,false,_sizeColDataY);

  // vector<int> valy = index(indexDataY,false,_sizeColDataY);
  bool equ = true;
  for(int j = 0; j < _systX[i].size(); j++)
  {
    // assert(_systX[i][j] <= _sizeColDataX);
    if(_cmi == false)
    {
      if((*_valX)(rowDataX,_systX[i][j]) != x[j])
      {
        equ = false;
      }
    }
    else
    {
      if((*_valX)(rowDataX,_systX[i][j]) != x[j])
      {
        equ = false;
      }
    }
  }

  for(int j=0;j<_systY[i].size();j++)
  {
    // assert(_systY[i][j]<= _sizeColDataY);
    if( valy[_systY[i][j]]!= y[j])
    {
      equ = false;
    }
  }

  if(equ)
  {
    return 1;
  }
  else
  {
    return -1;
  }
  delete[] x;
  delete[] y;
  delete[] valy;
}

int ITMatrix::getFeatureArraydeltaAlphYAlphX(int i, int indexX, int indexY, int indexDataX, int indexDataY)
{
  // assert(i<_systX.size());
  // assert(indexX < pow(_xAlphabet->rows(),_systX[i].size()));
  // assert(indexY < pow(_yAlphabet->rows(),_systY[i].size()));
  // assert(indexDataX < pow(_sizeX,_sizeColDataX));
  // assert(indexDataY < pow(_yAlphabet->rows(),_sizeColDataY));

  double* x    = new double[_systX[i].size()];
  double* y    = new double[_systY[i].size()];
  double* valy = new double[_sizeColDataY];
  double* valx = new double[_sizeColDataX];
  index(x, indexX, true,  _systX[i].size());
  index(y, indexY, false,  _systY[i].size());
  index(valy, indexDataY,false,_sizeColDataY);
  index(valx, indexDataX,true,_sizeColDataX);

  bool equ = true;

  for(int j = 0; j < _systX[i].size(); j++)
  {
    // assert(_systX[i][j] <= _sizeColDataX);
    if(valx[_systX[i][j]] != x[j])
    {
      equ = false;
    }
  }

  for(int j = 0; j < _systY[i].size(); j++)
  {
    assert(_systY[i][j]<= _sizeColDataY);
    if(valy[_systY[i][j]] != y[j])
    {
      equ = false;
    }
  }
  if(equ)
  {
    return 1;
  }
  else
  {
    return -1;
  }
  delete[] x;
  delete[] y;
  delete[] valy;
  delete[] valx;
}

// Berechnung der Alphabetreihe aus dem Index
// vector<int> ITMatrix::index(int index, bool x, int sizeCol)
void ITMatrix::index(double* array, int index, bool x, int sizeCol)
{
  int sizeAlph;
  // vector<int> row;
  double z;
  if(x)
  {
    // assert(index < pow(_xAlphabet->rows(),sizeCol));
    sizeAlph=_sizeX;
  }
  else
  {
    // assert(index < pow(_yAlphabet->rows(),sizeCol));
    sizeAlph=_sizeY;
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


