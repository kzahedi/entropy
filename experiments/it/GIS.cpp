#include "GIS.h"

#define EPSILON 0.00000001

// GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test, bool time, int seconds)
// :IT(eX, eY, aX, aY, lambdavalue, true)
GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
:IT(eX, eY, aX, aY,systX, systY, param, true)
{
  _param      = param;
  _expected   = new double**[_sizeSystX];
  for(int i=0; i<_sizeSystX; i++)
  {
      _expected[i]=new double*[(int) pow(_X->rows(),_systX[i].size())];
      for(int k=0; k<pow(_X->rows(),_systX[i].size()); k++)
      {
        _expected[i][k]= new double[(int) pow(_Y->rows(),_systY[i].size())];
        for(int l=0; l< pow(_Y->rows(),_systY[i].size());l++)
        {
          _expected[i][k][l]=0;
        }
      }
   }
  _normaliser = new double[_sizeSystX];
  _exponent= new double*[_sizeSystX];
  for(int i=0;i<_sizeSystX;i++){
	  _exponent[i]= new double[(int) pow(_Y->rows(),_sizeColValY)];
  }
  // cout << "Data X:" << endl << eX << endl << "Data Y: " << endl << eY << endl;
  if(param.time)
  {
    __gis(param.maxit, param.konv, param.test, param.seconds);
  }
  else
  {
    __gis(param.maxit, param.konv, param.test);
  }
}

GIS::~GIS()
{
  if(_observed!=NULL)
  {
    for(int i=0;i<_sizeSystX;i++)
    {
      for(int j=0;j<(int) pow(_X->rows(),_systX[i].size());j++)
      {
        delete [] _observed[i][j];
      }
     delete [] _observed[i];
  }
  delete [] _observed;
  }
  if(_expected != NULL){
     for(int i=0;i<_sizeSystX;i++){
         for(int k=0;k<(int) pow(_X->rows(),_systX[i].size());k++){
          delete [] _expected[i][k];
      }
      delete[] _expected[i];
     }
     delete[] _expected;
  }
  if(_exponent != NULL){
	  for(int i=0; i<_sizeSystX;i++){
		  delete [] _exponent[i];
	  }
    delete [] _exponent;
  }

  delete [] _normaliser;
}
//double GIS:: getconv(int i){
 // return _conv[i];
//}
//int GIS:: getsizeconv(){
//  return _conv.size();
//}


double GIS::__getFeatconst()
{
  double r = 0.0;
  for(int i=0; i< _sizeRowValX;i++) // i-th data row
  {
    for(int j=0; j< pow(_Y->rows(),_sizeColValY);j++) // y-alphabet
    {
      int v = (*_FM).getMatrixIndexFeat(i,j).size(); // the number of matching deltas
      if(v > r) r = v;
    }
  }
  return r;
}

void GIS::__getExpected()
{
	  for(int i=0; i<_sizeSystX; i++)
	  {
	    for(int k=0; k< (int) pow(_X->rows(),_systX[i].size()); k++)
	    {
	      for(int l=0; l< (int) pow(_Y->rows(),_systY[i].size());l++)
	      {
	        _expected[i][k][l]=0;
	      }
	    }
	  }
	  for(int xi=0; xi< _sizeRowValX; xi++)
	  {
	    for(int i=0;i<_sizeSystX;i++){
	      _normaliser[i]=0.0;
	    }
	    for(int yj=0; yj< pow(_Y->rows(),_sizeColValY); yj++)
	    {
	      for(int k=0; k< (*_FM).getMatrixIndexFeat(xi,yj).size();k++){
	    	int index=_FM->getMatrixIndexFeat(xi,yj)[k];
	    	_exponent[index][yj]=_FM->getFeatureArrayvalueAlphY(index,xi,yj);
	        _normaliser[index] += exp(_exponent[index][yj]);
	      }
	    }
	    for(int yj=0; yj< pow(_Y->rows(),_sizeColValY); yj++)
	    {
	      for(int k=0; k< (*_FM).getMatrixIndexFeat(xi,yj).size();k++)
	      {
	        int index=_FM->getMatrixIndexFeat(xi,yj)[k];
            _expected[index]
					 [(*_FM).getMatrixIndexdX(xi,yj)[k]]
	                 [(*_FM).getMatrixIndexdY(xi,yj)[k]] += exp(_exponent[index][yj])/_normaliser[index];
	      }
	    }
	   }
}
void GIS:: __gis(int maxit, double konv, bool test,int seconds){
    //constant c for delta
    double featconst = __getFeatconst();

    for(int k=0;k< _sizeSystX;k++){
      for(int i=0;i< pow(_Y->rows(),_sizeColValY);i++){
    	  _exponent[k][i]=0.0;
      }
    }
    _normaliser  = new double[_sizeSystX];
    _iterations  = 0;
    double l     = 1;
    double utime = 0;
    time_t befor;
    time_t after;
    while(utime<seconds )
    {
      befor=time(NULL);
      __calculateIteration(featconst, test);
      after=time(NULL);
      utime+= difftime(after,befor);
  //    cout << utime << endl;
    }
}

void GIS::__gis(int maxit, double konv, bool test)
{
  //constant c for delta
  double featconst = __getFeatconst();

  for(int k=0;k< _systX.size();k++){
    _exponent[k] = new double[(int) pow(_Y->rows(),_sizeColValY)]; // lambda_i * f_i
    for(int i=0;i< pow(_Y->rows(),_sizeColValY);i++){
  	  _exponent[k][i]=0.0;
    }
  }
  double l    = 1;
  _iterations = 0;
  while(_iterations < maxit && fabs(l) >= konv)
  {
    l = __calculateIteration(featconst, test);
  }
} 

double GIS::__calculateIteration(double featconst, bool test)
{
  double l = 0;
  __getExpected();
  for(int feat=0; feat < _sizeSystX; feat++)
  {
    // jedes delta hat ein x_i und ein y_j
    for(int deltai=0; deltai < pow(_X->rows(),_systX[feat].size()); deltai++)
    {
      for(int deltaj=0; deltaj < pow(_Y->rows(),_systY[feat].size()); deltaj++)
      {
        double oldl = (*_FM).getFeatureArraylambda(feat, deltai, deltaj);
        double newl = 0.0;
        if(fabs(_expected[feat][deltai][deltaj]) < EPSILON)
        {
          _expected[feat][deltai][deltaj] = 0.01; // TODO check if other values might be better
        }
        if(fabs(_observed[feat][deltai][deltaj]) > EPSILON)
        {
          // TODO 0.4 as learning rate parameter
          newl = oldl + 0.4*(1.0/featconst) *log(_observed[feat][deltai][deltaj]/_expected[feat][deltai][deltaj]);
        }
        else
        {
          newl = 0.0;
        }
        (*_FM).setFeatureArraylambda(feat,deltai,deltaj,newl);
        l+=fabs((_observed[feat][deltai][deltaj]-_expected[feat][deltai][deltaj]));
      }
    }
  }
  _iterations++;
  cout << _iterations << endl;
  cout << l << endl;
  if(test){
    _conv.push_back(l);
  }
  return l;

}

int GIS:: getIterations()
{
  return _iterations;
}
double GIS:: getconv(int i){
  return _conv[i];
}

int GIS:: getsizeconv(){
  return _conv.size();
}

