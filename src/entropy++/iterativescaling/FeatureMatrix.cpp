#include "FeatureMatrix.h"

#include <entropy++/powi.h>

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
  _sizeMatrixAlphabetY = powi(_yAlphabet->rows(),_sizeColDataY);
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
    for(int j=0;j<_sizeMatrixAlphabetY;j++)
    {
      for(int k=0;k<3 ;k++){
        _mat[i][j][k].clear();
      }
    }
  }
  delete _mat;
}

// ivector FeatureMatrix::getMatrixIndexFeat(int i, int j)
// {
  // ivector indX = _mat[i][j][0]; // TODO copying a vector can be expensive
  // return indX;
// }

// ivector FeatureMatrix::getMatrixIndexdX(int i, int j)
// {
  // ivector dindX = _mat[i][j][1]; // TODO copying a vector can be expensive
  // return dindX;
// }

// ivector FeatureMatrix::getMatrixIndexdY(int i, int j)
// {
  // ivector dindY = _mat[i][j][2]; // TODO copying a vector can be expensive
  // return dindY;
// }

int FeatureMatrix::getMatrixIndexFeatSize(int i, int j)
{
  return _mat[i][j][0].size();
}

int FeatureMatrix::getMatrixIndexFeatValue(int i, int j, int k)
{
  return _mat[i][j][0][k];
}
void FeatureMatrix::getMatrixIndexFeat(ivector &r, int i, int j)
{
  r = _mat[i][j][0]; // TODO copying a vector can be expensive
}

void FeatureMatrix::getMatrixIndexdX(ivector &r, int i, int j)
{
  r = _mat[i][j][1]; // TODO copying a vector can be expensive
}

int FeatureMatrix::getMatrixIndexdXValue(int i, int j, int k)
{
<<<<<<< Updated upstream
  assert(i < _sizeRowDataX && j < _sizeMatrixAlphabetY);
  ivector indX = _mat[i][j][0]; // TODO copying a vector can be expensive
  return indX;
=======
  return _mat[i][j][1][k]; // TODO copying a vector can be expensive
>>>>>>> Stashed changes
}

void FeatureMatrix::getMatrixIndexdY(ivector &r, int i, int j)
{
<<<<<<< Updated upstream
  assert(i < _sizeRowDataX  && j< _sizeMatrixAlphabetY);
  ivector dindX = _mat[i][j][1]; // TODO copying a vector can be expensive
  return dindX;
=======
  r = _mat[i][j][2]; // TODO copying a vector can be expensive
>>>>>>> Stashed changes
}

int FeatureMatrix::getMatrixIndexdYValue(int i, int j, int k)
{
<<<<<<< Updated upstream
  assert(i < _sizeRowDataX  && j < _sizeMatrixAlphabetY);
  ivector dindY = _mat[i][j][2]; // TODO copying a vector can be expensive
  return dindY;
=======
  return _mat[i][j][2][k]; // TODO copying a vector can be expensive
>>>>>>> Stashed changes
}

void FeatureMatrix::__getMatrix(double valuelambda)
{
  ivvector V(3,ivector(0));
  dvector x;
  dvector y;
  _mat = new ivvector*[_sizeRowDataX];
  for(int i = 0; i < _sizeRowDataX; i++)
  {
    _mat[i]= new ivvector[_sizeMatrixAlphabetY];
    for(int j = 0; j < _sizeMatrixAlphabetY; j++)
    {
      _mat[i][j] = V;
    }
  }

  for(int i = 0; i < _sizeRowDataX; i++)
  {
    for(int j = 0; j < _sizeMatrixAlphabetY; j++)
    {
      for(int feat = 0; feat < _systX.size(); feat++)
      {
        int DI = powi(_xAlphabet->rows(), _systX[feat].size());
        int DJ = powi(_yAlphabet->rows(), _systY[feat].size());
        cout << "  " << _sizeRowDataX * _sizeAlphY * _systX.size() * DI * DJ << endl;
        for(int deltai = 0; deltai < DI; deltai++)
        {
          for(int deltaj = 0; deltaj < DJ; deltaj++)
          {
            if(getDeltaAlphY(feat, deltai, deltaj, i, j) !=-1)
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
