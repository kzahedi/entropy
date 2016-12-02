#include "InstanceMatrix.h"

#include <entropy++/powi.h>

using namespace entropy::iterativescaling;

InstanceMatrix::InstanceMatrix(ULContainer *xData,
                               ULContainer *yData,
                               ULContainer *xAlphabet,
                               ULContainer *yAlphabet,
                               ivvector systX,
                               ivvector systY,
                               double valuelambda)
  : ITMatrix(xData,yData,xAlphabet,yAlphabet,systX, systY,valuelambda)
{
  __getMatrix(valuelambda);
}
InstanceMatrix::~InstanceMatrix()
{
  delete[] _featureArray;

  for(int i = 0; i < _systX.size(); i++)
  {
    int J = powi(_sizeX,_systX[i].size());
    int K = powi(_sizeY,_systY[i].size());
    for(int j = 0; j < J; j++)
    {
      for(int k = 0; k < K; k++)
      {
        _mat[i][j][k].clear();
      }
      delete[] _mat[i][j];
    }
    delete[] _mat[i];
  }
  delete[] _mat;
}

int InstanceMatrix::getInstanceMatrixX(int feat, int deltai, int deltaj, int k)
{
  assert(feat < _systX.size());
  assert(deltai < powi(_sizeX,_systX[feat].size()));
  assert(deltaj< powi(_sizeY,_systY[feat].size()));
  return _mat[feat][deltai][deltaj][0][k];
}

int InstanceMatrix::getInstanceMatrixY(int feat, int deltai, int deltaj, int k)
{
  assert(feat < _systX.size());
  assert(deltai < powi(_sizeX,_systX[feat].size()));
  assert(deltaj < powi(_sizeY,_systY[feat].size()));
  return _mat[feat][deltai][deltaj][1][k];
}
int InstanceMatrix::getInstanceMatrixSize(int feat, int deltai, int deltaj)
{
  assert(feat < _systX.size());
  assert(deltai < powi(_sizeX,_systX[feat].size()));
  assert(deltaj < powi(_sizeY,_systY[feat].size()));
  return _mat[feat][deltai][deltaj][0].size();

}
// Anzahl der deltas varriert
void InstanceMatrix::__getMatrix(double valuelambda)
{
  //Matrix erstellen
  ivvector V(2,ivector(0));
  _mat = new ivvector**[_systX.size()];
  for(int i = 0; i < _systX.size(); i++)
  {
    int K = (int)powi(_sizeX,_systX[i].size());
    int L = (int)powi(_sizeY,_systY[i].size());
    _mat[i] = new ivvector*[K];
    for(int k = 0; k < K; k++)
    {
      _mat[i][k]= new ivvector[L];
      for(int l = 0; l < L; l++)
      {
        _mat[i][k][l] = V;
      }
    }
  }

  //Matrix fuellen
  int Y = powi(_sizeY,_sizeColDataY);
  for(int feat = 0; feat < _systX.size(); feat++)
  {
    int DI = powi(_sizeX,_systX[feat].size());
    int DJ = powi(_sizeY,_systY[feat].size());
    for(int delti = 0; delti < DI; delti++)
    {
      for(int deltj = 0; deltj < DJ; deltj++)
      {
        for(int xi = 0; xi < _sizeRowDataX; xi++)
        {
          for(int y = 0; y < Y; y++)
          {
            if(getDeltaAlphY(feat,delti,deltj,xi,y)==1)
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
