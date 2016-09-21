#include "FeatureMatrix.h"

FeatureMatrix::FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double lambdavalue)
  :ITMatrix(eX,eY,aX,aY,lambdavalue)
{
  __getMatrix(lambdavalue);
}

FeatureMatrix::FeatureMatrix():ITMatrix()
{
  __getMatrix(0);
}

FeatureMatrix:: ~FeatureMatrix()
{
  for(int i=0; i<_sizeColValX;i++)
  {
    delete[] _FA[i];
  }
  delete _FA;

  for(int i=0;i<_sizeRowValX;i++)
  {
    for(int j=0;j<_sizeY;j++)
    {
      _mat[i][j].clear();
    }
  }
  delete _mat;
}

vector<int> FeatureMatrix::getMatrixIndexX(int i, int j)
{
  assert(i<_sizeRowValX && j< _sizeRowValY );
  vector<int> indX = _mat[i][j][0]; // TODO copying a vector can be expensive
  return indX;
}

vector<int> FeatureMatrix::getMatrixIndexY(int i, int j)
{
  assert(i<_sizeRowValX && j< _sizeY );
  vector<int> indY = _mat[i][j][1]; // TODO copying a vector can be expensive
  return indY;
}

vector<int> FeatureMatrix:: getMatrixIndexdX(int i,int j)
{
  assert(i<_sizeRowValX && j< _sizeY );
  vector<int> dindX = _mat[i][j][2]; // TODO copying a vector can be expensive
  return dindX;
}

vector<int> FeatureMatrix:: getMatrixIndexdY(int i,int j)
{
  assert(i<_sizeRowValX && j< _sizeY );
  vector<int> dindY = _mat[i][j][3]; // TODO copying a vector can be expensive
  return dindY;
}

void FeatureMatrix::__getMatrix(double valuelambda)
{
  vector<vector<int> > V(4,vector<int>(0));
  _mat = new vector<vector<int> >*[_sizeRowValX];
  for(int i=0;i<_sizeRowValX;i++)
  {
    _mat[i]= new vector<vector<int> >[_sizeY];
    for(int j=0;j<_sizeY;j++)
    {
      _mat[i][j]= V;
    }
  }
  for(int i=0;i<_sizeRowValX;i++)
  {
    for(int j=0;j<_sizeY;j++)
    {
      for(int varFeati=0;varFeati<_sizeColValX;varFeati++)
      {
        for(int varFeatj=0;varFeatj<_sizeColValY;varFeatj++)
        {
          if(_FA[varFeati][varFeatj].value((*_valX)(i,varFeati),(*_Y)(j,0))!=0)
          {
            for(int deltai=0; deltai<_sizeX; deltai++ )
            {
              for(int deltaj=0; deltaj<_sizeY; deltaj++)
              {
                if(_FA[varFeati][varFeatj].delta((*_X).get(deltai,0),
                                                 (*_Y).get(deltaj,0),
                                                 (*_valX).get(i, varFeati),
                                                 (*_Y).get(j,0)) != -1)
                {
                  _mat[i][j][0].push_back(varFeati);
                  _mat[i][j][1].push_back(varFeatj);
                  _mat[i][j][2].push_back(deltai);
                  _mat[i][j][3].push_back(deltaj);
                }
              }
            }
          }
        }
      }
    }
  }
}

