#include "FeatureMatrix.h"

using namespace entropy::iterativescaling;

FeatureMatrix::FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<vector<int> > systX, vector<vector<int> > systY,double lambdavalue)
  :ITMatrix(eX,eY,aX,aY,systX, systY,lambdavalue)
{
  _sizeAlphY = pow(_Y->rows(),_sizeColValY);
  __getMatrix(lambdavalue);
}
FeatureMatrix::FeatureMatrix(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY,double lambdavalue)
  :ITMatrix(eX,eY,aX,aY,systX,systY,lambdavalue)
{
  _sizeAlphY = pow(_Y->rows(),_sizeColValY);
  __getMatrix(lambdavalue);
}
FeatureMatrix::FeatureMatrix():ITMatrix()
{
  __getMatrix(0);
}

FeatureMatrix:: ~FeatureMatrix()
{
  delete _FA;

  for(int i=0;i<_sizeRowValX;i++)
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
vector<int> FeatureMatrix::getMatrixIndexFeat(int i,int j)
{
  assert(i<_sizeRowValX && j <_sizeAlphY);
  vector<int> indX = _mat[i][j][0]; // TODO copying a vector can be expensive
  return indX;
}

vector<int> FeatureMatrix:: getMatrixIndexdX(int i,int j)
{
  assert(i<_sizeRowValX  && j< _sizeAlphY );
  vector<int> dindX = _mat[i][j][1]; // TODO copying a vector can be expensive
  return dindX;
}

vector<int> FeatureMatrix:: getMatrixIndexdY(int i,int j)
{
  assert(i<_sizeRowValX  && j< _sizeAlphY );
  vector<int> dindY = _mat[i][j][2]; // TODO copying a vector can be expensive
  return dindY;
}

void FeatureMatrix::__getMatrix(double valuelambda)
{
  vector<vector<int> > V(3,vector<int>(0));
  vector<double > x;
  vector<double > y;
  _mat = new vector<vector<int> >*[_sizeRowValX];
  for(int i=0;i<_sizeRowValX;i++)
  {
    _mat[i]= new vector<vector<int> >[_sizeAlphY];
    for(int j=0;j<_sizeAlphY;j++)
    {
      _mat[i][j]= V;
    }
  }
  int j=3;
  for(int i=0;i<_sizeRowValX;i++)
  {
    for(int j=0;j<_sizeAlphY;j++)
    {
      for(int feat=0;feat<_systX.size();feat++)
      {
        for(int deltai=0; deltai<pow(_X->rows(),_systX[feat].size()); deltai++ )
        {
          for(int deltaj=0; deltaj<pow(_Y->rows(),_systY[feat].size()); deltaj++)
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

