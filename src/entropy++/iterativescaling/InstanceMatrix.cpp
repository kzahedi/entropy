#include "InstanceMatrix.h"

using namespace entropy::iterativescaling;

//die Alphabetwerte als double
InstanceMatrix:: InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,ivvector systX, ivvector systY, double valuelambda)
  :ITMatrix(eX,eY,aX,aY,systX, systY,valuelambda)
{
  _getMatrix(valuelambda);
}
//die Alphabetwerte als unsigned long
InstanceMatrix:: InstanceMatrix(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,ivvector systX, ivvector systY, double valuelambda)
  :ITMatrix(eX,eY,aX,aY,systX, systY,valuelambda)
{
  _getMatrix(valuelambda);
}

InstanceMatrix:: ~InstanceMatrix()
{
  delete [] _featureArray;

  for(int i=0;i<_systX.size();i++)
  {
    for(int j=0; j<  pow(_sizeX,_systX[i].size());j++)
    {
      for(int k=0; k< pow(_sizeY,_systY[i].size()); k++)
      {
        _mat[i][j][k].clear();
      }
      delete [] _mat[i][j];
    }
    delete [] _mat[i];
  }
  delete [] _mat;
}

ivector InstanceMatrix::getInstanceMatrixX(int feat, int deltai, int deltaj)
{
  assert(feat<_systX.size());
  assert(deltai<pow(_sizeX,_systX[feat].size()) && deltaj< pow(_sizeY,_systY[feat].size()));
  return _mat[feat][deltai][deltaj][0];
}
ivector InstanceMatrix::getInstanceMatrixY(int feat, int deltai, int deltaj)
{
  assert(feat<_systX.size());
  assert(deltai<pow(_sizeX,_systX[feat].size()) && deltaj< pow(_sizeY,_systY[feat].size()));
  return _mat[feat][deltai][deltaj][1];
}


// Anyahl der deltas varriert

void InstanceMatrix::_getMatrix(double valuelambda)
{
  //Matrix erstellen
  ivvector V(2,ivector(0));
  _mat = new ivvector**[_systX.size()];
  for(int i=0;i<_systX.size();i++)
  {
    _mat[i]= new ivvector*[(int) pow(_sizeX,_systX[i].size())];
    for(int k=0;k<pow(_sizeX,_systX[i].size());k++)
    {
      _mat[i][k]= new ivvector[(int) pow(_sizeY,_systY[i].size())];
      for(int l=0;l<pow(_sizeY,_systY[i].size());l++)
      {
        _mat[i][k][l]=V;
      }
    }
  }
  //Matrix fuellen
  for(int feat=0;feat<_systX.size();feat++)
  {
    for(int delti=0; delti<pow(_sizeX,_systX[feat].size());delti++)
    {
      for(int deltj=0; deltj <pow(_sizeY,_systY[feat].size()); deltj++)
      {
        for(int xi=0; xi< _sizeRowDataX; xi++)
        {
          for(int y=0; y<pow(_sizeY,_sizeColDataY); y++)
          {
            if(getFeatureArraydeltaAlphY(feat,delti,deltj,xi,y)==1)
            {
              _mat[feat][delti][deltj][0].push_back(xi);
              _mat[feat][delti][deltj][1].push_back(y);
            }
          }
        }
      }
    }
  }
}
