#include "GIS.h"

#define EPSILON 0.00000001

GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test, bool time, int seconds)
:IT(eX, eY, aX, aY, lambdavalue, true)
{
  _exponent   = new double[_sizeY];
  _normaliser = 0.0;
  _expected   = new double***[_sizeColValX];
  for(int i=0; i<_sizeColValX; i++)
  {
    _expected[i]=new double**[_sizeColValY];
    for( int j=0;j< _sizeColValY;j++)
    {
      _expected[i][j]=new double*[_sizeX];
      for(int k=0; k< _sizeX; k++)
      {
        _expected[i][j][k]= new double[_sizeY];
        for(int l=0; l< _sizeY;l++)
        {
          _expected[i][j][k][l]=0;
        }
      }
    }
  }

  // cout << "Data X:" << endl << eX << endl << "Data Y: " << endl << eY << endl;

  if(time)
  {
    __gis(maxit, konv, test, seconds);
  }
  else
  {
    __gis(maxit, konv, test);
  }
}

GIS::~GIS()
{
  if(_observed!=NULL)
  {
    for(int i=0;i<_sizeColValX;i++)
    {
      for(int j=0;j<_sizeColValY;j++)
      {
        for(int k=0;k<_sizeX;k++)
        {
           delete [] _observed[i][j][k];
      }
      delete [] _observed[i][j];
     }
     delete [] _observed[i];
  }
  delete [] _observed;
  }
  if(_expected != NULL){
     for(int i=0;i<_sizeColValX;i++){
        for(int j=0;j<_sizeColValY;j++){
         for(int k=0;k<_sizeX;k++){
          delete [] _expected[i][j][k];
       }
       delete [] _expected[i][j];
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
      int v = (*_FM).getMatrixIndexX(i,j).size(); // the number of matching deltas
      if(v > r) r = v;
    }
  }
  return r;
}

void GIS::__getExpected()
{
  for(int i=0; i<_sizeColValX; i++)
  {
    for( int j=0;j< _sizeColValY;j++)
    {
      for(int k=0; k< _sizeX; k++)
      {
        for(int l=0; l< _sizeY;l++)
        {
          _expected[i][j][k][l]=0;
        }
      }
    }
  }
  for(int Feati=0; Feati< _sizeColValX; Feati++)
  {
    for(int Featj=0; Featj< _sizeColValY; Featj++)
    {
      for(int xi=0; xi< _sizeRowValX; xi++)
      {
        _normaliser = 0.0;
        for(int yj=0; yj< _sizeY; yj++)
        {
          _exponent[yj]=_FM->getFeatureArrayvalue(Feati,Featj,(*_valX)(xi,Feati), (*_Y)(yj,0));
          _normaliser += exp(_exponent[yj]);
        }
        for(int yj=0; yj< _sizeY; yj++)
        {
          for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++)
          {
            if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj)
            {
              if((*_FM).getFeatureArraydelta(Feati,
                                             Featj,
                                             (*_FM).getMatrixIndexdX(xi,yj)[k],
                                             (*_FM).getMatrixIndexdY(xi,yj)[k],
                                             (*_valX)(xi,Feati),
                                             (*_Y)(yj,0)) == 1)
              {
              _expected[Feati]
                       [Featj]
                       [(*_FM).getMatrixIndexdX(xi,yj)[k]]
                       [(*_FM).getMatrixIndexdY(xi,yj)[k]] += exp(_exponent[yj])/_normaliser;
              }
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
  for(int Feati=0; Feati < _sizeColValX; Feati++)
  {
    for(int Featj=0; Featj < _sizeColValY; Featj++)
    {
      // jedes delta hat ein x_i und ein y_j
      for(int deltai=0; deltai < _sizeX; deltai++)
      {
        for(int deltaj=0; deltaj < _sizeY; deltaj++)
        {
          double oldl = (*_FM).getFeatureArraylambda(Feati, Featj, deltai, deltaj);
          double newl = 0.0;
          if(fabs(_expected[Feati][Featj][deltai][deltaj]) < EPSILON)
          {
            _expected[Feati][Featj][deltai][deltaj] = 0.01; // TODO check if other values might be better
          }
          if(fabs(_observed[Feati][Featj][deltai][deltaj]) > EPSILON)
          {
            // TODO 0.1 as learning rate parameter
            newl = oldl + 0.1 * (1.0/featconst) *
              log(_observed[Feati][Featj][deltai][deltaj]/_expected[Feati][Featj][deltai][deltaj]);
          }
          else
          {
            newl = 0.0;
          }
          (*_FM).setFeatureArraylambda(Feati,Featj,deltai,deltaj,newl);
          l+=fabs((_observed[Feati][Featj][deltai][deltaj]-_expected[Feati][Featj][deltai][deltaj]));
        }
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

