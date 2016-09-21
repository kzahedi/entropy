#include "IT.h"

IT::IT(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue, bool GIS)
{
  assert(eX.rows() == eY.rows());
  _gis         = GIS;
  _valX        = &eX;
  _valY        = &eY;
  _X           = &aX;
  _Y           = &aY;
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  _sizeColValY = _valY->columns();
  _sizeColValX = _valX->columns();
  _sizeRowValX = _valX->rows();
  _sizeRowValY = _valY->rows();

  if(_gis) // gis and csgis require different feature matrices
  {
    _FM = new FeatureMatrix(*_valX,*_valY,*_X,*_Y,lambdavalue);
  }
  else
  {
    _IM = new InstanceMatrix(*_valX,*_valY,*_X,*_Y,lambdavalue);
  }
  _observed=__getobs();
}

IT::IT(int ColValY, DContainer &eX, DContainer &aX, DContainer &aY)
{
  _gis         = true;
  _X           = &aX;
  _Y           = &aY;
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  _valX        = &eX;
  _sizeColValY = ColValY;
  _sizeColValX = (*_valX).columns();
  _sizeRowValX = (*_valX).rows();
  _valY        = new DContainer(_sizeRowValX,ColValY);
  _sizeRowValY = 0;
  _FM          = new FeatureMatrix(*_valX,*_valY,*_X,*_Y,1);
  _IM          = NULL;
  _observed    = NULL;
}

// returns p(y_j = valY | x_i = valX)
double IT::prop(int Feati, int Featj, double ValX, double ValY)
{
  assert(Feati<_sizeColValX && Featj<_sizeColValY);
  double norm     = 0;
  double exponent = 0;
  if(_gis)
  {
    exponent = exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,ValY));
  }
  else
  {
    exponent = exp((*_IM).getFeatureArrayvalue(Feati,Featj,ValX,ValY));
  }

  for(int yi = 0; yi < _sizeY; yi++)
  {
    if(_gis)
    {
      norm += exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
    }
    else
    {
      norm += exp((*_IM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
    }
  }
  return exponent/norm;
}

// p(y | x)
// double IT::prop(int rowX, vector<vector<double> > Y, int rowY)
double IT::prop(int rowX, vector<vector<double> >& Y, int rowY)
{
  double feat     = 0;
  double featnorm = 0;
  double norm     = 0;
  double exponent = 0;
  for(int Featx=0; Featx< _sizeColValX; Featx++)
  {
    for(int Featy=0; Featy< _sizeColValY; Featy++)
    {
      if(_gis)
      {
        feat += (*_FM).getFeatureArrayvalue(Featx,Featy,(*_valX)(rowX,Featx), Y[rowY][Featy]);
      }
      else
      {
        feat += (*_IM).getFeatureArrayvalue(Featx,Featy,(*_valX)(rowX,Featx), Y[rowY][Featy]); // Y = _alphY
      }
    }
  }
  exponent = exp(feat);
  for(int yi=0;yi<Y.size();yi++)
  {
    for(int Featx=0;Featx< _sizeColValX;Featx++)
    {
      for(int Featy=0;Featy< _sizeColValY;Featy++)
      {
        if(_gis)
        {
          featnorm += (*_FM).getFeatureArrayvalue(Featx,Featy,(*_valX)(rowX,Featx),Y[yi][Featy]);
        }
        else
        {
          featnorm += (*_IM).getFeatureArrayvalue(Featx,Featy,(*_valX)(rowX,Featx),Y[yi][Featy]);
        }
      }
    }

    norm     += exp(featnorm);
    featnorm  = 0;
  }
  return exponent/norm;
}

// P(x)
double IT::propm(vector<vector<double> > X, int rowX, vector<vector<double> >& Y){
  double z=0;
  double feat=0;
  for(int y=0;y<Y.size();y++)
  {
    for(int Featx=0;Featx< _sizeColValX;Featx++)
    {
      for(int Featy=0;Featy< _sizeColValY;Featy++)
      {
        if(_gis)
        {
          feat+=(*_FM).getFeatureArrayvalue(Featx,Featy,X[rowX][Featx],Y[y][Featy]);
        }
        else
        {
          feat+=(*_IM).getFeatureArrayvalue(Featx,Featy,X[rowX][Featx],Y[y][Featy]);
        }
      }
    }
    z+=exp(feat);
    feat=0;
  }
  double n;
  double nexp=0;
  for(int x=0;x<X.size();x++)
  {
    for(int y=0;y<Y.size();y++)
    {
      for(int Featx=0;Featx< _sizeColValX;Featx++)
      {
        for(int Featy=0;Featy< _sizeColValY;Featy++)
        {
          if(_gis)
          {
            feat+=(*_FM).getFeatureArrayvalue(Featx,Featy,X[x][Featx],Y[y][Featy]);
          }
          else
          {
            feat+=(*_IM).getFeatureArrayvalue(Featx,Featy,X[x][Featx],Y[y][Featy]);
          }
        }
      }
      nexp +=exp(feat);
      feat=0;
    }
    n+=nexp;
    nexp=0;
  }
  return z/n;
}

double IT::getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY)
{
  assert(Feati<_sizeColValX && Featj<_sizeColValY);
  double lambda = 0.0;
  if(_gis)
  {
    lambda = _FM->getFeatureArraylambda(Feati,Featj, ilambdaX, ilambdaY);
  }
  else
  {
    lambda = _IM->getFeatureArraylambda(Feati,Featj, ilambdaX, ilambdaY);
  }
  return lambda;
}

// get observed
double**** IT::__getobs()
{
  _observed = new double***[_sizeColValX];
  for(int i=0; i<_sizeColValX; i++)
  {
    _observed[i]=new double**[_sizeColValY];
    for( int j=0;j< _sizeColValY;j++)
    {
      _observed[i][j]=new double*[_sizeX];
      for(int k=0; k< _sizeX; k++)
      {
        _observed[i][j][k]= new double[_sizeY];
        for(int l=0; l< _sizeY;l++)
        {
          _observed[i][j][k][l] = 0.0;
        }
      }
    }
  }

  //vector observed
  for(int i=0;i<_sizeRowValX;i++ )
  {
    for(int feati=0; feati< _sizeColValX;feati++)
    {
      for(int featj=0;featj< _sizeColValY; featj++)
      {
        for(int delti = 0; delti < _sizeX; delti++)
        {
          for(int deltj=0; deltj<_sizeY; deltj++)
          {
            if(_gis)
            {
              if(_FM->getFeatureArraydelta(feati,featj,delti,deltj,(*_valX)(i,feati),(*_valY)(i,featj))==1)
              {
                _observed[feati][featj][delti][deltj]++;
              }
            }
            else
            {
              if(_IM->getFeatureArraydelta(feati,featj,delti,deltj,(*_valX)(i,feati),(*_valY)(i,featj))==1)
              {
                _observed[feati][featj][delti][deltj]++;
              }
            }
          }
        }
      }
    }
  }
  return _observed;
}

void IT::setFeatureArraylambda(int feati, int featj, int ilambdaX, int ilambdaY,double valuelambda){
  assert(feati<_sizeColValX && featj<_sizeColValY);
  assert(_gis == true);
  _FM->setFeatureArraylambda(feati,featj,ilambdaX,ilambdaY,valuelambda);
}


