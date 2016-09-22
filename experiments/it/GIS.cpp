#include "GIS.h"

#define EPSILON 0.00000001

// GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test, bool time, int seconds)
// :IT(eX, eY, aX, aY, lambdavalue, true)
GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
:IT(eX, eY, aX, aY,systX, systY, param, true)
{
  _param      = param;
  _exponent   = new double[_sizeY];
  _normaliser = 0.0;
  _expected   = new double**[_systX.size()];
  for(int i=0; i<_sizeColValX; i++)
  {
      _expected[i]=new double*[_sizeX];
      for(int k=0; k< _sizeX; k++)
      {
        _expected[i][k]= new double[_sizeY];
        for(int l=0; l< _sizeY;l++)
        {
          _expected[i][k][l]=0;
        }
      }
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
    for(int i=0;i<_systX.size();i++)
    {
      for(int j=0;j<_sizeX;j++)
      {
        delete [] _observed[i][j];
      }
     delete [] _observed[i];
  }
  delete [] _observed;
  }
  if(_expected != NULL){
     for(int i=0;i<_systX.size();i++){
         for(int k=0;k<_sizeX;k++){
          delete [] _expected[i][k];
      }
      delete[] _expected[i];
     }
     delete[] _expected;
  }
  if(_exponent != NULL){
    delete [] _exponent;
  }
  _conv.clear();
}
double GIS:: getconv(int i){
  return _conv[i];
}
int GIS:: getsizeconv(){
  return _conv.size();
}


double GIS::__getFeatconst()
{
  double r = 0.0;
  for(int i=0; i< _sizeRowValX;i++) // i-th data row 
  {
    for(int j=0; j< _sizeY;j++) // y-alphabet
    {
      int v = (*_FM).getMatrixIndexFeat(i,j).size(); // the number of matching deltas
      if(v > r) r = v;
    }
  }
  return r;
}

void GIS::__getExpected()
{
  for(int i=0; i<_systX.size(); i++)
  {
    for(int k=0; k< _sizeX; k++)
    {
      for(int l=0; l< _sizeY;l++)
      {
        _expected[i][k][l]=0;
      }
    }

  }
  for(int Feati=0; Feati< _systX.size(); Feati++)
  {
    for(int xi=0; xi< _sizeRowValX; xi++)
    {
      _normaliser = 0.0;
      for(int yj=0; yj< pow(_Y->rows(),_sizeColValY); yj++)
      {
        _exponent[yj]=_FM->getFeatureArrayvalueAlphY(Feati,xi,yj);
          _normaliser += exp(_exponent[yj]);
      }
      for(int yj=0; yj< pow(_Y->rows(),_sizeColValY); yj++)
      {
        for(int k=0; k< (*_FM).getMatrixIndexFeat(xi,yj).size();k++)
        {
          if((*_FM).getMatrixIndexFeat(xi,yj)[k]==Feati)
          {
            if((*_FM).getFeatureArraydeltaAlphY(Feati, (*_FM).getMatrixIndexdX(xi,yj)[k], (*_FM).getMatrixIndexdY(xi,yj)[k], xi, yj) == 1)
            {
              _expected[Feati]
				       [(*_FM).getMatrixIndexdX(xi,yj)[k]]
                       [(*_FM).getMatrixIndexdY(xi,yj)[k]] += exp(_exponent[yj])/_normaliser;
            }
          }
        }
      }
    }
  }
}
void GIS:: __gis(int maxit, double konv, bool test,int seconds){

    //constant c for delta
    double featconst = __getFeatconst();

    for(int k=0;k< _sizeY;k++){
      _exponent[k] = 0.0; // lambda_i * f_i
    }
    _normaliser  = 0.0;
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
    }
}

void GIS::__gis(int maxit, double konv, bool test)
{
  //constant c for delta
  double featconst = __getFeatconst();

  for(int k=0; k < _sizeY; k++)
  {
    _exponent[k]=0.0; // lambda_i * f_i
  }
  _normaliser = 0.0;
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
  for(int Feati=0; Feati < _systX.size(); Feati++)
  {
    // jedes delta hat ein x_i und ein y_j
    for(int deltai=0; deltai < _sizeX; deltai++)
    {
      for(int deltaj=0; deltaj < _sizeY; deltaj++)
      {
        double oldl = (*_FM).getFeatureArraylambda(Feati, deltai, deltaj);
        double newl = 0.0;
        if(fabs(_expected[Feati][deltai][deltaj]) < EPSILON)
        {
          _expected[Feati][deltai][deltaj] = 0.01; // TODO check if other values might be better
        }
        if(fabs(_observed[Feati][deltai][deltaj]) > EPSILON)
        {
          // TODO 0.1 as learning rate parameter
          newl = oldl + 0.1 * (1.0/featconst) *log(_observed[Feati][deltai][deltaj]/_expected[Feati][deltai][deltaj]);
        }
        else
        {
          newl = 0.0;
        }
        (*_FM).setFeatureArraylambda(Feati,deltai,deltaj,newl);
        l+=fabs((_observed[Feati][deltai][deltaj]-_expected[Feati][deltai][deltaj]));
      }
    }
  }
  _iterations++;
  if(test){
    _conv.push_back(l);
  }
  return l;
}

int GIS:: getIterations()
{
  return _iterations;
}

