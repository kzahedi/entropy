#include "InstanceMatrix.h"

InstanceMatrix:: InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, double valuelambda)
  :ITMatrix(eX,eY,aX,aY,systX, systY,valuelambda)
{
  _getMatrix(valuelambda);
}

InstanceMatrix:: ~InstanceMatrix()
{
  delete _FA;

  for(int i=0;i<_systX.size();i++)
  {
     _mat[i].clear();
  }
  delete _mat;
}

vector<int> InstanceMatrix::getInstanceMatrixDeltaX(int feat){
  assert(feat<_systX.size());
  return _mat[feat][0];
}
vector<int> InstanceMatrix::getInstanceMatrixDeltaY(int feat){
  assert(feat<_systX.size());
  return _mat[feat][1];
}
vector<int> InstanceMatrix::getInstanceMatrixX(int feat)
{
  assert(feat<_systX.size());
  return _mat[feat][2];
}
vector<int> InstanceMatrix::getInstanceMatrixY(int feat)
{
  assert(feat<_systX.size());
  return _mat[feat][3];
}


// Anyahl der deltas varriert

void InstanceMatrix::_getMatrix(double valuelambda)
{
  vector<vector<int> > V(4,vector<int>(0));
  _mat = new vector<vector<int> >[_systX.size()];
  for(int i=0;i<_systX.size();i++)
  {
        _mat[i]=V;
  }

  for(int feat=0;feat<_systX.size();feat++)
  {
    for(int delti=0; delti<pow(_X->rows(),_systX[feat].size());delti++)
    {
      for(int deltj=0; deltj < pow(_Y->rows(),_systY[feat].size()); deltj++)
      {
        for(int xi=0; xi< _sizeRowValX; xi++)
        {
          for(int y=0; y< pow(_Y->rows(),_sizeColValY); y++)
            {
              if(getFeatureArraydeltaAlphY(feat, delti, deltj, xi, y)==1)
              {
                _mat[feat][0].push_back(delti);
                _mat[feat][1].push_back(deltj);
                _mat[feat][2].push_back(xi);
                _mat[feat][3].push_back(y);
              }
            }
          }
        }
      }

  }
}
