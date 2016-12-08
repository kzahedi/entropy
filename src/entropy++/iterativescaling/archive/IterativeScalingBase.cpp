#include "IterativeScalingBase.h"

#include <entropy++/powi.h>

using namespace entropy::iterativescaling;

IterativeScalingBase::IterativeScalingBase(ULContainer *xData,
                                           ULContainer *yData,
                                           ULContainer *xAlphabet,
                                           ULContainer *yAlphabet,
                                           ivvector systX,
                                           ivvector systY,
                                           IsParameter param,
                                           bool useFeatures)
{
  assert(xData->rows() == yData->rows());
  _param            = param;
  _xData            = NULL;
  _yData            = NULL;
  ULContainer *valX = xData;
  ULContainer *valY = yData;
  _xAlphabet        = xAlphabet;
  _yAlphabet        = yAlphabet;
  _sizeSystX        = systX.size();
  _sizeX            = _xAlphabet->rows();
  _sizeY            = _yAlphabet->rows();
  _sizeColDataY     = valY->columns();
  _sizeColDataX     = valX->columns();
  _sizeRowDataX     = valX->rows();
  _sizeRowDataY     = valY->rows();
  _systX            = systX;
  _systY            = systY;

  if(useFeatures) // isGis and csisGis require different feature matrices
  {
    _imatrix = new FeatureMatrix(valX,
                                 valY,
                                 _xAlphabet,
                                 _yAlphabet,
                                 systX,
                                 systY,
                                 param.lambdavalue);
  }
  else
  {
    _imatrix = new InstanceMatrix(valX,
                                  valY,
                                  _xAlphabet,
                                  _yAlphabet,
                                  systX,
                                  systY,
                                  param.lambdavalue);
  }
  _observed = __getobs();
}

IterativeScalingBase::IterativeScalingBase(int ColDataY,
                                           ULContainer *xData,
                                           ULContainer *xAlphabet,
                                           ULContainer *yAlphabet,
                                           ivvector   systX,
                                           ivvector   systY)
{
  _systX        = systX;
  _systY        = systY;
  _xAlphabet    = xAlphabet;
  _yAlphabet    = yAlphabet;
  _sizeSystX    = systX.size();
  _sizeX        = _xAlphabet->rows();
  _sizeY        = _yAlphabet->rows();
  _xData        = xData;
  _sizeColDataY = ColDataY;
  _sizeColDataX = _xData->columns();
  _sizeRowDataX = _xData->rows();
  _yData        = new ULContainer(_sizeRowDataX,ColDataY);
  _sizeRowDataY = _sizeRowDataX;
  _imatrix      = new FeatureMatrix(_xData,
                                    _yData,
                                    _xAlphabet,
                                    _yAlphabet,
                                    systX,
                                    systY,
                                    1);
  _observed     = NULL;
}

// p(y | x)
// double IterativeScalingBase::prop(int rowX,  int indexY)
double IterativeScalingBase::prop(int rowX, int indexY)
{
  double featexp  = 0;
  double featnorm = 0;
  double norm     = 0;
  double exponent = 1;

  for(int feat = 0; feat < _systX.size(); feat++)
  {
    featexp += _imatrix->getValueAlphY(feat,rowX, indexY);
  }

  int YI = powi(_yAlphabet->rows(),_sizeColDataY);
  for(int yi = 0; yi < YI; yi++)
  {
    for(int feat = 0; feat < _systX.size(); feat++)
    {
      featnorm += _imatrix->getValueAlphY(feat, rowX, yi);
    }
    norm     += exp(featnorm-featexp);
    featnorm  = 0;
  }
  return exponent/norm;
}

// p(y | x)
// double IterativeScalingBase::prop(int indexX, int indexY)
// double IterativeScalingBase::propAlphX(int indexX, int rowY)
double IterativeScalingBase::p_x_c_y(int indexX, int rowY)
{
  double featexp  = 0;
  double featnorm = 0;
  double norm     = 0;
  double exponent = 1;
  for(int feat = 0; feat < _systX.size(); feat++)
  {
    featexp += _imatrix->getValueAlphYAlphX(feat, indexX, rowY);
  }
  int YI = powi(_yAlphabet->rows(),_sizeColDataY);
  for(int yi = 0; yi < YI; yi++)
  {
    for(int feat = 0; feat < _systX.size(); feat++)
    {
      featnorm += _imatrix->getValueAlphYAlphX(feat, indexX, yi);
    }
    norm     += exp(featnorm-featexp);
    featnorm  = 0;
  }
  return exponent/norm;
}

// P(x)
double IterativeScalingBase::propm(int rowX)
{
  double max      = 0;
  double z        = 0;
  double curr;
  double exponent = 0;
  int YI          = powi(_yAlphabet->rows(),_sizeColDataY);
  for(int y = 0; y < YI; y++)
  {
    for(int feat = 0; feat < _systX.size(); feat++)
    {
      curr = _imatrix->getValueAlphYAlphX(feat,rowX,y);
      if(fabs(curr) > fabs(max)) max = curr; 
    }
  }
  for(int y = 0; y < YI; y++)
  {
    for(int feat = 0; feat < _systX.size(); feat++)
    {
      exponent += _imatrix->getValueAlphYAlphX(feat,rowX,y);
    }
    z        += exp(exponent - max);  // = 0 fuer grosse lambda
    exponent  = 0;
  }
  double n;
  double nexp = 0;
  int    XI   = powi(_xAlphabet->rows(),_sizeColDataX);
  for(int x = 0; x < XI; x++)
  {
    for(int y = 0; y < YI; y++)
    {
      for(int feat = 0; feat < _systX.size(); feat++)
      {
        exponent += _imatrix->getValueAlphYAlphX(feat,x,y);
      }
      nexp     += exp(exponent-max);
      exponent  = 0;
    }
    n    += nexp;
    nexp  = 0;
  }
  return z/n;
}

double IterativeScalingBase::getFeatureArraylambda(int Feati, int ilambdaX, int ilambdaY)
{
  // kontrolle der Parameter ilambdaX und ilambdaY findet in Feature statt
  assert(Feati<_systX.size());
  return _imatrix->getLambda(Feati, ilambdaX, ilambdaY);
}

// get observed
double*** IterativeScalingBase::__getobs()
{
  //erzeugen von observed
  _observed = new double**[_sizeSystX];
  for(int i = 0; i < _systX.size(); i++)
  {
    int K = (int)powi(_xAlphabet->rows(), _systX[i].size());
    int L = (int)powi(_yAlphabet->rows(), _systY[i].size());
    _observed[i] = new double*[K];
    for(int k = 0; k < K; k++)
    {
      _observed[i][k] = new double[L];
      for(int l = 0; l < L; l++)
      {
        _observed[i][k][l] = 0.0;
      }
    }
  }

  // fuellen von observed
  for(int i = 0; i < _sizeRowDataX; i++)
  {
    for(int feat = 0; feat< _sizeSystX; feat++)
    {
      int DI = powi(_xAlphabet->rows(), _systX[feat].size());
      int DJ = powi(_yAlphabet->rows(), _systY[feat].size());
      for(int delti = 0; delti < DI; delti++)
      {
        for(int deltj=0; deltj < DJ; deltj++)
        {
          if(_imatrix->getDelta(feat, delti, deltj, i, i) == 1)
          {
            _observed[feat][delti][deltj]++;
          }
        }
      }
    }
  }
  return _observed;
}

IterativeScalingBase::~IterativeScalingBase()
{
  if(_observed!=NULL)
  {
    for(int i = 0; i <_sizeSystX; i++)
    {
      int K = powi(_xAlphabet->rows(), _systX[i].size());
      for(int k = 0; k < K; k++)
      {
        delete[] _observed[i][k];
      }
      delete[] _observed[i];
    }
    delete[] _observed;
  }
  for(int i = 0; i < _systX.size(); i++) _systX[i].clear();
  _systX.clear();

  for(int j = 0; j < _systY.size(); j++) _systY[j].clear();
  _systY.clear();
}

//veraendern der lambdas im Algorithmus funktioniert mit einer Methode in ITMatrix, hier nur zum erzeugen von Testwerten
void IterativeScalingBase::setFeatureArraylambda(int feati, int ilambdaX, int ilambdaY,double valuelambda)
{
  assert(feati<_systX.size());
  _imatrix->setLambda(feati, ilambdaX, ilambdaY, valuelambda);
}

// Berechnung der Alphabetreihe aus dem Index
dvector IterativeScalingBase::index(int index, bool x, int sizeCol)
{
  int sizeAlph;
  dvector row;
  double z;
  if(x)
  {
    assert(index < powi(_xAlphabet->rows(), sizeCol));
    sizeAlph = _sizeX;
  }
  else
  {
    assert(index < powi(_yAlphabet->rows(), sizeCol));
    sizeAlph = _sizeY;
  }
  for(int i = sizeCol; i>0 ; i--)
  {
    z = index % ((int)powi(sizeAlph, i));
    z = z / (powi(sizeAlph,i-1));
    int j = (int)z;
    if(x)
    {
      row.push_back((*_xAlphabet)(j, 0));
    }
    else
    {
      row.push_back((*_yAlphabet)(j, 0));
    }
  }
  return row;
}
