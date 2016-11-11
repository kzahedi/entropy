#include "IterativeScaling.h"

using namespace entropy::iterativescaling;

IterativeScaling::IterativeScaling(DContainer &xData,
                                   DContainer &yData,
                                   DContainer &xAlphabet,
                                   DContainer &yAlphabet,
                                   vector<vector<int> > systX,
                                   vector<vector<int> > systY,
                                   IsParameter param,
                                   bool isGis)
{
  assert(xData.rows() == yData.rows());
  _param       = param;
  _isGis       = isGis;
  _xData       = &xData;
  _yData       = &yData;
  _xAlphabet   = &xAlphabet;
  _yAlphabet   = &yAlphabet;
  _sizeSystX   = systX.size();
  _sizeX       = _xAlphabet->rows();
  _sizeY       = _yAlphabet->rows();
  _sizeColDataY = _yData->columns();
  _sizeColDataX = _xData->columns();
  _sizeRowDataX = _xData->rows();
  _sizeRowDataY = _yData->rows();
  _systX       = systX;
  _systY       = systY;

  if(_isGis) // isGis and csisGis require different feature matrices
  {
    _FM = new FeatureMatrix(*_xData, *_yData, *_xAlphabet, *_yAlphabet, systX, systY, param.lambdavalue);
  }
  else
  {
    _IM = new InstanceMatrix(*_xData, *_yData, *_xAlphabet, *_yAlphabet, systX, systY, param.lambdavalue);
  }
  _observed = __getobs();
}
// Konstruktor fuer Alphabete mit unsigned long
IterativeScaling::IterativeScaling(ULContainer &xData,
                                   ULContainer &yData,
                                   DContainer &xAlphabet,
                                   DContainer &yAlphabet,
                                   vector<vector<int> > systX,
                                   vector<vector<int> > systY,
                                   IsParameter param,
                                   bool isGis)
{
  assert(xData.rows() == yData.rows());
  _param            = param;
  _isGis            = isGis;
  _xData            = NULL;
  _yData            = NULL;
  ULContainer *valX = &xData;
  ULContainer *valY = &yData;
  _xAlphabet        = &xAlphabet;
  _yAlphabet        = &yAlphabet;
  _sizeSystX        = systX.size();
  _sizeX            = _xAlphabet->rows();
  _sizeY            = _yAlphabet->rows();
  _sizeColDataY     = valY->columns();
  _sizeColDataX     = valX->columns();
  _sizeRowDataX     = valX->rows();
  _sizeRowDataY     = valY->rows();
  _systX            = systX;
  _systY            = systY;

  if(_isGis) // isGis and csisGis require different feature matrices
  {
    _FM = new FeatureMatrix(*valX,*valY,*_xAlphabet,*_yAlphabet, systX, systY, param.lambdavalue);
  }
  else
  {
    _IM = new InstanceMatrix(*valX,*valY,*_xAlphabet,*_yAlphabet, systX, systY, param.lambdavalue);
  }

  _observed = __getobs();
}

IterativeScaling::IterativeScaling(int ColDataY,
                                   DContainer &xData,
                                   DContainer &xAlphabet,
                                   DContainer &yAlphabet,
                                   vector<vector<int> > systX,
                                   vector<vector<int> > systY)
{
  _isGis        = true;
  _systX        = systX;
  _systY        = systY;
  _xAlphabet    = &xAlphabet;
  _yAlphabet    = &yAlphabet;
  _sizeSystX    = systX.size();
  _sizeX        = _xAlphabet->rows();
  _sizeY        = _yAlphabet->rows();
  _xData        = &xData;
  _sizeColDataY = ColDataY;
  _sizeColDataX = (*_xData).columns();
  _sizeRowDataX = (*_xData).rows();
  _yData        = new DContainer(_sizeRowDataX,ColDataY);
  _sizeRowDataY = _sizeRowDataX;
  _FM           = new FeatureMatrix(*_xData,*_yData,*_xAlphabet,*_yAlphabet,systX,systY,1);
  _IM           = NULL;
  _observed     = NULL;
}

// p(y | x)
// double IterativeScaling::prop(int rowX,  int indexY)
double IterativeScaling::prop(int rowX, int indexY)
{
  double featexp  = 0;
  double featnorm = 0;
  double norm     = 0;
  double exponent = 1;

  for(int feat = 0; feat < _systX.size(); feat++)
  {
    if(_isGis)
    {
      featexp += (*_FM).getFeatureArrayvalueAlphY(feat,rowX, indexY);
    }
    else
    {
      featexp += (*_IM).getFeatureArrayvalueAlphY(feat,rowX, indexY);
    }
  }

  for(int yi=0;yi<pow(_yAlphabet->rows(),_sizeColDataY);yi++)
  {
    for(int feat=0; feat< _systX.size(); feat++)
    {
      if(_isGis)
      {
        featnorm += (*_FM).getFeatureArrayvalueAlphY(feat,rowX,yi);
      }
      else
      {
        featnorm += (*_IM).getFeatureArrayvalueAlphY(feat,rowX,yi);
      }
    }
    norm     += exp(featnorm-featexp);
    featnorm  = 0;
  }

  return exponent/norm;
}

// p(y | x)
// double IterativeScaling::prop(int indexX, int indexY)
double IterativeScaling::propAlphX(int indexX, int rowY)
{
  double featexp  = 0;
  double featnorm = 0;
  double norm     = 0;
  double exponent = 1;
  for(int feat=0; feat< _systX.size(); feat++)
  {
    if(_isGis)
    {
      featexp += (*_FM).getFeatureArrayvalueAlphYAlphX(feat, indexX, rowY);
    }
    else
    {
      featexp += (*_IM).getFeatureArrayvalueAlphYAlphX(feat, indexX, rowY);
    }
  }
  for(int yi=0;yi<pow(_yAlphabet->rows(),_sizeColDataY);yi++)
  {
    for(int feat=0; feat< _systX.size(); feat++)
    {
      if(_isGis)
      {
        featnorm += (*_FM).getFeatureArrayvalueAlphYAlphX(feat,indexX,yi);
      }
      else
      {
        featnorm += (*_IM).getFeatureArrayvalueAlphYAlphX(feat,indexX,yi);
      }
    }
    norm     += exp(featnorm-featexp);
    featnorm  = 0;
  }
  return exponent/norm;
}

// P(x)
double IterativeScaling::propm(int rowX)
{
  double z        = 0.0;
  double exponent = 0.0;

  for(int y=0;y<pow(_yAlphabet->rows(),_sizeColDataY);y++)
  {
    for(int feat=0;feat<  _systX.size();feat++)
    {
      if(_isGis)
      {
        exponent+=(*_FM).getFeatureArrayvalueAlphYAlphX(feat,rowX,y);
      }
      else
      {
        exponent+=(*_IM).getFeatureArrayvalueAlphYAlphX(feat,rowX,y);
      }

    }
    z+=exp(exponent);
    exponent=0;
  }

  double n    = 0.0;
  double nexp = 0.0;

  for(int x=0; x <pow(_xAlphabet->rows(),_sizeColDataX);x++)
  {
    for(int y=0;y<pow(_yAlphabet->rows(),_sizeColDataY);y++)
    {
      for(int feat=0;feat< _systX.size();feat++)
      {
        if(_isGis)
        {
          exponent+=(*_FM).getFeatureArrayvalueAlphYAlphX(feat,x,y);
        }
        else
        {
          exponent+=(*_IM).getFeatureArrayvalueAlphYAlphX(feat,x,y);
        }
      }
      nexp +=exp(exponent);
      exponent=0;
    }
    n+=nexp;
    nexp=0;
  }
  return z/n;
}

double IterativeScaling::getFeatureArraylambda(int Feati, int ilambdaX, int ilambdaY)
{
  assert(Feati<_systX.size()); //kontrolle der Parameter ilambdaX und ilambdaY findet in Feature statt
  double lambda = 0.0;
  if(_isGis)
  {
    lambda = _FM->getFeatureArraylambda(Feati, ilambdaX, ilambdaY);
  }
  else
  {
    lambda = _IM->getFeatureArraylambda(Feati, ilambdaX, ilambdaY);
  }
  return lambda;
}

// get observed
double*** IterativeScaling::__getobs()
{
  //erzeugen von observed
  _observed = new double**[_sizeSystX];
  for(int i=0; i<_systX.size(); i++)
  {
    _observed[i]=new double*[(int) pow(_xAlphabet->rows(),_systX[i].size())];
    for(int k=0; k< pow(_xAlphabet->rows(),_systX[i].size()); k++)
    {
      _observed[i][k]= new double[(int) pow(_yAlphabet->rows(),_systY[i].size())];
      for(int l=0; l< pow(_yAlphabet->rows(),_systY[i].size());l++)
      {
        _observed[i][k][l] = 0.0;
      }
    }
  }
  // fuellen von observed
  for(int i=0;i<_sizeRowDataX;i++ )
  {
    for(int feat=0; feat< _sizeSystX;feat++)
    {
      for(int delti = 0; delti<pow(_xAlphabet->rows(),_systX[feat].size()); delti++)
      {
        for(int deltj=0; deltj<pow(_yAlphabet->rows(),_systY[feat].size()); deltj++)
        {
          if(_isGis)
          {
            if(_FM->getFeatureArraydelta(feat,delti,deltj,i,i)==1)
            {
              _observed[feat][delti][deltj]++;
            }
          }
          else
          {
            if(_IM->getFeatureArraydelta(feat,delti,deltj,i,i)==1)
            {
              _observed[feat][delti][deltj]++;
            }
          }
        }
      }
    }
  }
  return _observed;
}
IterativeScaling::~IterativeScaling(){
  if(_observed!=NULL)
  {
    for(int i=0;i<_sizeSystX;i++)
    {
      for(int k=0;k<pow(_xAlphabet->rows(),_systX[i].size());k++)
      {
        delete [] _observed[i][k];
      }
      delete [] _observed[i];
    }
    delete [] _observed;
  }
  for(int i=0; i< _systX.size();i++)
  {
    _systX[i].clear();
  }
  _systX.clear();
  for(int j=0; j<_systY.size();j++)
  {
    _systY[j].clear();
  }
  _systY.clear();
}
//veraendern der lambdas im Algorithmus funktioniert mit einer Methode in ITMatrix, hier nur zum erzeugen von Testwerten
void IterativeScaling::setFeatureArraylambda(int feati, int ilambdaX, int ilambdaY,double valuelambda){
  assert(feati<_systX.size());
  assert(_isGis == true);
  _FM->setFeatureArraylambda(feati,ilambdaX,ilambdaY,valuelambda);
}

// Berechnung der Alphabetreihe aus dem Index
vector<double> IterativeScaling::index(int index,bool x,int sizeCol)
{
  int sizeAlph;
  vector<double> zeile;
  double z;
  if(x){
    assert(index<pow(_xAlphabet->rows(),sizeCol));
    sizeAlph=_sizeX;
  }
  else{
    assert(index<pow(_yAlphabet->rows(),sizeCol));
    sizeAlph=_sizeY;
  }
  for(int i = sizeCol; i>0 ; i--)
  {
    z = index % (int) (pow(sizeAlph,i));
    z = z/ (pow(sizeAlph,i-1));
    int j= (int) z;
    if(x){
      zeile.push_back((*_xAlphabet)(j,0));
    }
    else{
      zeile.push_back((*_yAlphabet)(j,0));
    }
  }
  return zeile;
}
