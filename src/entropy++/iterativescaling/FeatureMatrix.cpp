#include "FeatureMatrix.h"

using namespace entropy::iterativescaling;

FeatureMatrix::FeatureMatrix(ULContainer *xData,
                             ULContainer *yData,
                             ULContainer *xAlphabet,
                             ULContainer *yAlphabet,
                             ivvector systX,
                             ivvector systY,
                             double lambdavalue)
  : ITMatrix(xData,yData,xAlphabet,yAlphabet,systX,systY,lambdavalue)
{
  _sizeAlphY = pow(_yAlphabet->rows(),_sizeColDataY);
  __getMatrix(lambdavalue);
}

FeatureMatrix::FeatureMatrix()
  : ITMatrix()
{
  __getMatrix(0);
}

FeatureMatrix::~FeatureMatrix()
{
  delete _featureArray;

  for(int i=0;i<_sizeRowDataX;i++)
  {
    for(int j=0;j<_sizeAlphY;j++)
    {
      for(int k=0;k<3 ;k++){
        _mat[i][j][k].clear();
      }
    }
  }
  delete _mat;
}

ivector FeatureMatrix::getMatrixIndexFeat(int i, int j)
{
  assert(i < _sizeRowDataX && j < _sizeAlphY);
  ivector indX = _mat[i][j][0]; // TODO copying a vector can be expensive
  return indX;
}

ivector FeatureMatrix::getMatrixIndexdX(int i, int j)
{
  assert(i < _sizeRowDataX  && j< _sizeAlphY);
  ivector dindX = _mat[i][j][1]; // TODO copying a vector can be expensive
  return dindX;
}

ivector FeatureMatrix::getMatrixIndexdY(int i, int j)
{
  assert(i < _sizeRowDataX  && j < _sizeAlphY);
  ivector dindY = _mat[i][j][2]; // TODO copying a vector can be expensive
  return dindY;
}

void FeatureMatrix::__getMatrix(double valuelambda)
{
  ivvector V(3,ivector(0));
  dvector x;
  dvector y;
  _mat = new ivvector*[_sizeRowDataX];
  for(int i = 0; i < _sizeRowDataX; i++)
  {
    _mat[i]= new ivvector[_sizeAlphY];
    for(int j = 0; j < _sizeAlphY; j++)
    {
      _mat[i][j] = V;
    }
  }

  for(int i = 0; i < _sizeRowDataX; i++)
  {
    for(int j = 0; j < _sizeAlphY; j++)
    {
      for(int feat = 0; feat < _systX.size(); feat++)
      {
        int DI = pow(_xAlphabet->rows(), _systX[feat].size());
        int DJ = pow(_yAlphabet->rows(), _systY[feat].size());
        for(int deltai = 0; deltai < DI; deltai++)
        {
          for(int deltaj = 0; deltaj < DJ; deltaj++)
          {
            if(getFeatureArraydeltaAlphY(feat, deltai, deltaj, i, j) !=-1)
            {
              _mat[i][j][0].push_back(feat);
              _mat[i][j][1].push_back(deltai);
              _mat[i][j][2].push_back(deltaj);
            }
          }
        }
      }
    }
  }
}
