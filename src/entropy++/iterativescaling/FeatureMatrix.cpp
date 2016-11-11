#include "FeatureMatrix.h"

using namespace entropy::iterativescaling;

FeatureMatrix::FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, ivvector systX, ivvector systY,double lambdavalue)
  :ITMatrix(eX,eY,aX,aY,systX, systY,lambdavalue)
{
  _sizeAlphY = pow(_yAlphabet->rows(),_sizeColDataY);
  __getMatrix(lambdavalue);
}

FeatureMatrix::FeatureMatrix(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,ivvector systX, ivvector systY,double lambdavalue)
  :ITMatrix(eX,eY,aX,aY,systX,systY,lambdavalue)
{
  _sizeAlphY = pow(_yAlphabet->rows(),_sizeColDataY);
  __getMatrix(lambdavalue);
}
FeatureMatrix::FeatureMatrix():ITMatrix()
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

ivector FeatureMatrix::getMatrixIndexFeat(int i,int j)
{
  assert(i<_sizeRowDataX && j <_sizeAlphY);
  ivector indX = _mat[i][j][0]; // TODO copying a vector can be expensive
  return indX;
}

ivector FeatureMatrix::getMatrixIndexdX(int i,int j)
{
  assert(i<_sizeRowDataX  && j< _sizeAlphY );
  ivector dindX = _mat[i][j][1]; // TODO copying a vector can be expensive
  return dindX;
}

ivector FeatureMatrix::getMatrixIndexdY(int i,int j)
{
  assert(i<_sizeRowDataX  && j< _sizeAlphY );
  ivector dindY = _mat[i][j][2]; // TODO copying a vector can be expensive
  return dindY;
}

void FeatureMatrix::__getMatrix(double valuelambda)
{
  ivvector V(3,ivector(0));
  vector<double > x;
  vector<double > y;
  _mat = new ivvector*[_sizeRowDataX];
  for(int i=0;i<_sizeRowDataX;i++)
  {
    _mat[i]= new ivvector[_sizeAlphY];
    for(int j=0;j<_sizeAlphY;j++)
    {
      _mat[i][j]= V;
    }
  }
  int j=3;
  for(int i=0;i<_sizeRowDataX;i++)
  {
    for(int j=0;j<_sizeAlphY;j++)
    {
      for(int feat=0;feat<_systX.size();feat++)
      {
        for(int deltai=0; deltai<pow(_xAlphabet->rows(),_systX[feat].size()); deltai++ )
        {
          for(int deltaj=0; deltaj<pow(_yAlphabet->rows(),_systY[feat].size()); deltaj++)
          {
            if(getFeatureArraydeltaAlphY(feat,deltai,deltaj,i,j) !=-1)
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
