#include "InstanceMatrix.h"

InstanceMatrix:: InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, double valuelambda)
  :ITMatrix(eX,eY,aX,aY,systX, systY,valuelambda)
{
  _getMatrix(valuelambda);
}

InstanceMatrix:: InstanceMatrix(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, double valuelambda)
  :ITMatrix(eX,eY,aX,aY,systX, systY,valuelambda)
{
	cout << " vor _getMatrix " << endl;
  _getMatrix(valuelambda);
}

InstanceMatrix:: ~InstanceMatrix()
{
  delete [] _FA;

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

vector<int> InstanceMatrix::getInstanceMatrixX(int feat, int deltai, int deltaj)
{
   assert(feat<_systX.size());
   assert(deltai<pow(_sizeX,_systX[feat].size()) && deltaj< pow(_sizeY,_systY[feat].size()));
	return _mat[feat][deltai][deltaj][0];
}
vector<int> InstanceMatrix::getInstanceMatrixY(int feat, int deltai, int deltaj)
{
  assert(feat<_systX.size());
  assert(deltai<pow(_sizeX,_systX[feat].size()) && deltaj< pow(_sizeY,_systY[feat].size()));
  return _mat[feat][deltai][deltaj][1];
}


// Anyahl der deltas varriert

void InstanceMatrix::_getMatrix(double valuelambda)
{
	cout << " instanceMatrix " << endl;
	cout << " hier " << endl;
	vector<vector<int> > V(2,vector<int>(0));
	_mat = new vector<vector<int> >**[_systX.size()];
	for(int i=0;i<_systX.size();i++)
	{
	  _mat[i]= new vector<vector<int> >*[(int) pow(_sizeX,_systX[i].size())];
	  for(int k=0;k<pow(_sizeX,_systX[i].size());k++)
	  {
	    _mat[i][k]= new vector<vector<int> >[(int) pow(_sizeY,_systY[i].size())];
	    for(int l=0;l<pow(_sizeY,_systY[i].size());l++)
	    {
	   	  _mat[i][k][l]=V;
        }
	  }
	}
	cout << " nach dem ersten Block " << endl;
	for(int feat=0;feat<_systX.size();feat++)
	{
		cout << " hier 1 " << endl;
      for(int delti=0; delti<pow(_sizeX,_systX[feat].size());delti++)
      {
    		cout << " hier 2 " << endl;
	    for(int deltj=0; deltj <pow(_sizeY,_systY[feat].size()); deltj++)
	    {
	    	cout << " hier 3 " << endl;
		  for(int xi=0; xi< _sizeRowValX; xi++)
		  {
				cout << " hier 4 " << endl;
		    for(int y=0; y<pow(_sizeY,_sizeColValY); y++)
		    {
		    	cout << " hier 5 " << endl;
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
