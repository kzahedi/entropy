#include "ITMatrix.h"

using namespace entropy::iterativescaling;

ITMatrix::ITMatrix(DContainer &xData,
                   DContainer &yData,
                   DContainer &xAlphabet,
                   DContainer &yAlphabet,
                   vector<vector<int> > systX,
                   vector<vector<int> > systY,
                   double lambdavalue)
{
  assert(systX.size()==systY.size());
  _valX         = &xData;
  _valY         = &yData;
  _xAlphabet    = &xAlphabet;
  _yAlphabet    = &yAlphabet;
  _systX        = systX;
  _systY        = systY;
  _sizeColDataY = _valY->columns();
  _sizeColDataX = _valX->columns();
  _sizeRowDataX = _valX->rows();
  _sizeRowDataY = _valY->rows();
  _sizeX        = _xAlphabet->rows();
  _sizeY        = _yAlphabet->rows();
  _cmi          = false;
  __featureArray(lambdavalue);
}

//die Alphabetwerte als unsigned long
ITMatrix::ITMatrix(ULContainer &xData,
                   ULContainer &yData,
                   DContainer &xAlphabet,
                   DContainer &yAlphabet,
                   vector<vector<int> > systX,
                   vector<vector<int> > systY,
                   double lambdavalue)
{
  assert(systX.size() == systY.size());
  _valX         = NULL;
  _valY         = NULL;
  _valXUL       = &xData;
  _valYUL       = &yData;
  _xAlphabet    = &xAlphabet;
  _yAlphabet    = &yAlphabet;
  _cmi          = true;
  _systX        = systX;
  _systY        = systY;
  _sizeColDataY = _valYUL->columns();
  _sizeColDataX = _valXUL->columns();
  _sizeRowDataX = _valXUL->rows();
  _sizeRowDataY = _valYUL->rows();
  _sizeX        = _xAlphabet->rows();
  _sizeY        = _yAlphabet->rows();
  __featureArray(lambdavalue);
}

ITMatrix::ITMatrix()
{
  _valX         = new DContainer(0,0);
  _valY         = new DContainer(0,0);
  _xAlphabet    = new DContainer(0,0);
  _yAlphabet    = new DContainer(0,0);
  _sizeColDataY = 0;
  _sizeColDataX = 0;
  _sizeRowDataX = 0;
  _sizeRowDataY = 0;
  _sizeX        = _xAlphabet->rows();
  _sizeY        = _yAlphabet->rows();
  __featureArray(0);
}

ITMatrix::~ITMatrix()
{
  delete _featureArray;
  for(int i=0;i<_systX.size();i++){
    _systX[i].clear();
  }
  _systX.clear();
  for(int j=0;j<_systY.size();j++){
    _systY[j].clear();
  }
  _systY.clear();
}

double  ITMatrix::getFeatureArraylambda(int i,int ilambdaX, int ilambdaY)
{
  assert(i<_systX.size());  //ilambdaX und ilambdaY werden in Feature geprueft
  double lambda=_featureArray[i].getLambda(ilambdaX,ilambdaY);
  return lambda;
}

double ITMatrix::getFeatureArrayvalue(int feat,int rowX,int rowY)
{
  assert(feat<_systX.size());
  double val=0.0;
  for(int i=0; i<  pow(_xAlphabet->rows(),_systX[feat].size());i++){
    for(int j=0;j <pow(_yAlphabet->rows(),_systY[feat].size());j++){
      val+= getFeatureArraydelta(feat,i,j,rowX,rowY)*_featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}
double ITMatrix::getFeatureArrayvalueAlphY(int feat,int rowX,int indexY){
  assert(feat<_systX.size());
  double val=0.0;
  for(int i=0; i<  pow(_xAlphabet->rows(),_systX[feat].size());i++){
    for(int j=0;j <pow(_yAlphabet->rows(),_systY[feat].size());j++){
      val+= getFeatureArraydeltaAlphY(feat,i,j,rowX,indexY)*_featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}
double ITMatrix::getFeatureArrayvalueAlphYAlphX(int feat,int indexX,int indexY){
  assert(feat<_systX.size());
  double val=0.0;
  for(int i=0; i<  pow(_xAlphabet->rows(),_systX[feat].size());i++){
    for(int j=0;j <pow(_yAlphabet->rows(),_systY[feat].size());j++){
      val+= getFeatureArraydeltaAlphYAlphX(feat,i,j,indexX,indexY)*_featureArray[feat].getLambda(i,j);
    }
  }
  return val;
}
void  ITMatrix::setFeatureArraylambda(int i, int ilambdaX, int ilambdaY,double valuelambda)
{
  assert(i<_systX.size());
  _featureArray[i].setLambda(ilambdaX,ilambdaY,valuelambda);
}
int   ITMatrix::getFeatureArraydelta(int i,int indexX, int indexY,int rowDataX, int rowDataY)
{
  assert(i<_systX.size());
  assert(indexX < pow(_xAlphabet->rows(),_systX[i].size()) && indexY<pow(_yAlphabet->rows(),_systY[i].size()));
  assert(rowDataX < _sizeRowDataX);
  assert(rowDataY < _sizeRowDataY);
  vector<double> x = index(indexX,true,_systX[i].size());
  vector<double> y = index(indexY,false,_systY[i].size());
  bool equ = true;
  for(int j=0; j<_systX[i].size();j++)
  {
    assert(_systX[i][j]<= _sizeColDataX);
    if(_cmi==false){
      if((*_valX)(rowDataX,_systX[i][j]) != x[j]){
        equ = false;
      }
    }
    else{
      if((*_valXUL)(rowDataX,_systX[i][j]) != x[j]){
        equ = false;
      }
    }
  }
  for(int j=0;j<_systY[i].size();j++){
    assert(_systY[i][j]<= _sizeColDataY);
    if(_cmi==false){
      if((*_valY)(rowDataY,_systY[i][j]) != y[j]){
        equ = false;
      }
    }
    else{
      if((*_valYUL)(rowDataY,_systY[i][j]) != y[j]){
        equ = false;
      }
    }

  }
  if(equ){
    return 1;
  }
  else{
    return -1;
  }

}
// Feature i, indizes der zu vergleichenden Werte aus den tatsaechlich auftetenden X und den moeglichen Y
int   ITMatrix::getFeatureArraydeltaAlphY(int i,int indexX, int indexY,int rowDataX, int indexDataY)
{
  assert(i<_systX.size());
  assert(indexX < pow(_xAlphabet->rows(),_systX[i].size()) && indexY<pow(_yAlphabet->rows(),_systY[i].size()));
  assert(rowDataX < _sizeRowDataX);
  assert(indexDataY < pow(_yAlphabet->rows(),_sizeColDataY));
  vector<double> x = index(indexX,true,_systX[i].size());
  vector<double> y = index(indexY,false,_systY[i].size());
  vector<double> valy = index(indexDataY,false,_sizeColDataY);
  bool equ = true;
  for(int j=0; j<_systX[i].size();j++)
  {
    assert(_systX[i][j]<= _sizeColDataX);
    if(_cmi==false){
      if((*_valX)(rowDataX,_systX[i][j]) != x[j]){
        equ = false;
      }
    }
    else{
      if((*_valXUL)(rowDataX,_systX[i][j]) != x[j]){
        equ = false;
      }
    }
  }
  for(int j=0;j<_systY[i].size();j++){
    assert(_systY[i][j]<= _sizeColDataY);
    if( valy[_systY[i][j]]!= y[j]){
      equ = false;
    }
  }
  if(equ){
    return 1;
  }
  else{
    return -1;
  }
}
// Feature i, indizes der zu vergleichenden Werte aus den moeglichen X und Y
int ITMatrix::getFeatureArraydeltaAlphYAlphX(int i,int indexX, int indexY,int indexDataX, int indexDataY)
{
  assert(i<_systX.size());
  assert(indexX < pow(_xAlphabet->rows(),_systX[i].size()) && indexY<pow(_yAlphabet->rows(),_systY[i].size()));
  assert(indexDataX < pow(_sizeX,_sizeColDataX));
  assert(indexDataY < pow(_yAlphabet->rows(),_sizeColDataY));
  vector<double> x = index(indexX,true, _systX[i].size());
  vector<double> y = index(indexY,false,_systY[i].size());
  vector<double> valy = index(indexDataY,false, _sizeColDataY);
  vector<double> valx = index(indexDataX,true,_sizeColDataX);
  bool equ = true;
  for(int j=0; j<_systX[i].size();j++)
  {
    assert(_systX[i][j]<= _sizeColDataX);
    if(valx[_systX[i][j]] != x[j]){
      equ = false;
    }
  }
  for(int j=0;j<_systY[i].size();j++){
    assert(_systY[i][j]<= _sizeColDataY);
    if( valy[_systY[i][j]]!= y[j]){
      equ = false;
    }
  }
  if(equ){
    return 1;
  }
  else{
    return -1;
  }
}
// Berechnung der Alphabetreihe aus dem Index
vector<double> ITMatrix::index(int index,bool x,int sizeCol)
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

//ein Array mit den benoetigten Features fuellen
void ITMatrix::__featureArray(double valuelambda)
{
  int sizeX = 0;
  int sizeY = 0;

  _featureArray = new Feature[_systX.size()];

  for(int m=0; m < _systX.size(); m++)
  {
    // Feature *K = new Feature(*_xAlphabet,*_yAlphabet,_systX[m].size(),_systY[m].size(),valuelambda);
    sizeX = pow(_xAlphabet->rows(), _systX[m].size());
    sizeY = pow(_yAlphabet->rows(), _systY[m].size());
    Feature *K = new Feature(sizeX, sizeY, valuelambda);
    _featureArray[m]= *K;
  }
}
