#include "FeatureMatrix.h"
//DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<vector<int> > systX, vector<vector<int> > systY, double lambdavalue
FeatureMatrix::FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<vector<int> > systX, vector<vector<int> > systY,double lambdavalue)
  :ITMatrix(eX,eY,aX,aY,systX, systY,lambdavalue)
{
	 cout << "hier 0002 " << endl;
  _sizeAlphY = pow(_Y->rows(),_sizeColValY);
  cout << "hier 0002 " << endl;
  __getMatrix(lambdavalue);
  cout << "hier 62 " << endl;

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
    for(int j=0;j<_sizeY;j++)
    {
      _mat[i][j].clear();
    }
  }
  delete _mat;
}
vector<int> FeatureMatrix::getMatrixIndexFeat(int i,int j)
{
  assert(i<_systX.size());
  vector<int> indX = _mat[i][j][0]; // TODO copying a vector can be expensive
  return indX;
}

vector<int> FeatureMatrix:: getMatrixIndexdX(int i,int j)
{
  assert(i<_sizeRowValX && j< _sizeY );
  vector<int> dindX = _mat[i][j][1]; // TODO copying a vector can be expensive
  return dindX;
}

vector<int> FeatureMatrix:: getMatrixIndexdY(int i,int j)
{
  assert(i<_sizeRowValX && j< _sizeY );
  vector<int> dindY = _mat[i][j][2]; // TODO copying a vector can be expensive
  return dindY;
}

void FeatureMatrix::__getMatrix(double valuelambda)
{
	 cout << "hier 0002 " << endl;
  vector<vector<int> > V(3,vector<int>(0));
  vector<double > x;
  vector<double > y;
  _mat = new vector<vector<int> >*[_sizeRowValX];
  for(int i=0;i<_sizeRowValX;i++)
  {
    _mat[i]= new vector<vector<int> >[_sizeAlphY];
    for(int j=0;j<_sizeY;j++)
    {
      _mat[i][j]= V;
    }
  }
  cout << _sizeAlphY << "                       hier " << endl;
  cout << "hier 1 " << endl;
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
              { cout << "if " << i << j << feat << deltai << deltaj <<  endl;
                _mat[i][j][0].push_back(feat);
                _mat[i][j][1].push_back(deltai);
                _mat[i][j][2].push_back(deltaj);
              }
              else{
            	  cout << "else " << i << j << feat << deltai << deltaj << endl;
              }
            }
         }
      }
    }
    cout << "hier fast " << endl;
  }
  cout << "hier ende " << endl;
}

