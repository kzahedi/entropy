#include "ITMatrix.h"

ITMatrix::ITMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<vector<int> > systX, vector<vector<int> > systY, double lambdavalue)
{
  _valX        = &eX;
  _valY        = &eY;
  _X           = &aX;
  _Y           = &aY;
  assert(systX.size()==systY.size());
  _systX       = systX;
  _systY       = systY;
  _sizeColValY = _valY->columns();
  _sizeColValX = _valX->columns();
  _sizeRowValX = _valX->rows();
  _sizeRowValY = _valY->rows();
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  FeatureArray(lambdavalue);
}

ITMatrix::ITMatrix()
{
  _valX        = new DContainer(0,0);
  _valY        = new DContainer(0,0);
  _X           = new DContainer(0,0);
  _Y           = new DContainer(0,0);
  _sizeColValY = 0;
  _sizeColValX = 0;
  _sizeRowValX = 0;
  _sizeRowValY = 0;
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  FeatureArray(0);
}

ITMatrix::~ITMatrix()
{
  delete _FA;
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
  assert(i<_systX.size());
  double lambda=_FA[i].getlambda(ilambdaX,ilambdaY);
  return lambda;
}

double ITMatrix::getFeatureArrayvalue(int feat,int rowX,int rowY)
{
 assert(feat<_systX.size());
 double val;
 for(int i=0; i<  pow(_X->rows(),_systX[feat].size());i++){
	 for(int j=0;j <pow(_Y->rows(),_systY[feat].size());j++){
		 val+= getFeatureArraydelta(feat,i,j,rowX,rowY)*_FA[feat].getlambda(i,j);
	 }
 }
 return val;
}
double ITMatrix::getFeatureArrayvalueAlphY(int feat,int rowX,int indexY){
	 assert(feat<_systX.size());
	 double val;
	 for(int i=0; i<  pow(_X->rows(),_systX[feat].size());i++){
		 for(int j=0;j <pow(_Y->rows(),_systY[feat].size());j++){
			 val+= getFeatureArraydeltaAlphY(feat,i,j,rowX,indexY)*_FA[feat].getlambda(i,j);
		     if(indexY == 0 && feat == 0 && rowX == 0){
		//    	 cout << " lambda " << _FA[feat].getlambda(i,j)<< "  " <<  getFeatureArraydeltaAlphY(feat,i,j,rowX,indexY) << " " << i << " " << j<< endl;
		     }
		 }
	 }
	 return val;
	 cout << " neu " << endl;
}
double ITMatrix::getFeatureArrayvalueAlphYAlphX(int feat,int indexX,int indexY){
	 assert(feat<_systX.size());
	 double val;
	 for(int i=0; i<  pow(_X->rows(),_systX[feat].size());i++){
		 for(int j=0;j <pow(_Y->rows(),_systY[feat].size());j++){
			 val+= getFeatureArraydeltaAlphYAlphX(feat,i,j,indexX,indexY);
		 }
	 }
	 return val;
}
void  ITMatrix::setFeatureArraylambda(int i, int ilambdaX, int ilambdaY,double valuelambda)
{
  assert(i<_systX.size());
  _FA[i].setlambda(ilambdaX,ilambdaY,valuelambda);
}
int   ITMatrix::getFeatureArraydelta(int i,int indexX, int indexY,int rowValX, int rowValY)
{
  assert(i<_systX.size());
  vector<double> x = index(indexX,true,_systX[i].size());
  vector<double> y = index(indexY,false,_systY[i].size());
  bool equ = true;
  for(int j=0; j<_systX[i].size();j++)
  {
	  if((*_valX)(rowValX,_systX[i][j]) != x[j]){
		  equ = false;
	  }
	//  cout << (*_valX)(rowValX,_systX[i][j]) << " " << _systX[i][j] << " " << i << " " <<  j << " " << x[_systX[i][j]] << "row " << rowValX << " " << rowValY << " " << equ <<  endl;
  }
  for(int j=0;j<_systY[i].size();j++){
	  if((*_valY)(rowValY,_systY[i][j]) != y[j]){
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
int   ITMatrix::getFeatureArraydeltaAlphY(int i,int indexX, int indexY,int rowValX, int indexValY)
{
  assert(i<_systX.size());
  //assert(idelta < _sizeX && jdelta < _sizeY);
  vector<double> x = index(indexX,true,_systX[i].size());
  vector<double> y = index(indexY,false,_systY[i].size());
  vector<double> valy = index(indexValY,false, _sizeColValY);
  bool equ = true;
  for(int j=0; j<_systX[i].size();j++)
  {
	  if((*_valX)(rowValX,_systX[i][j]) != x[j]){
		  equ = false;
	  }
	//  cout << (*_valX)(rowValX,_systX[i][j]) << " " << _systX[i][j] << " " << i << " " <<  j << " " << x[j] << "row " << rowValX << " " << indexValY << " " << equ <<  endl;
  }
  for(int j=0;j<_systY[i].size();j++){
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
// Feature i, indizes der zu vergleichenden Werte
int ITMatrix::getFeatureArraydeltaAlphYAlphX(int i,int indexX, int indexY,int indexValX, int indexValY)
{
  assert(i<_systX.size());
  //assert(idelta < _sizeX && jdelta < _sizeY);
  vector<double> x = index(indexX,true, _systX[i].size());
  vector<double> y = index(indexY,false,_systY[i].size());
  vector<double> valy = index(indexValY,false, _sizeColValY);
  vector<double> valx = index(indexValX,true,_sizeColValX);
  bool equ = true;
  for(int j=0; j<_systX[i].size();j++)
  {
	  if(valx[_systX[i][j]] != x[j]){
		  equ = false;
	  }
  }
  for(int j=0;j<_systY[i].size();j++){
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
	  sizeAlph=_sizeX;
  }
  else{
	  sizeAlph=_sizeY;
  }
  for(int i = sizeCol; i>0 ; i--)
  {
	  z = index % (int) (pow(sizeAlph,i));
	  z = z/ (pow(sizeAlph,i-1));
	  int j= (int) z;
	  if(x){
		  zeile.push_back((*_X)(j,0));
	  }
	  else{
		  zeile.push_back((*_Y)(j,0));
	  }
  }
  return zeile;
}
void ITMatrix::FeatureArray( double valuelambda)
{
  _FA = new Feature[_systX.size()];
  for(int m=0; m< _systX.size(); m++)
  { //(DContainer &aX, DContainer &aY,int colValX, int colValY, int systX, int systXsize,int systYsize , double valuelambda);
      Feature *K = new Feature(*_X,*_Y,_sizeColValX,_sizeColValY,m,_systX[m].size(),_systY[m].size(),valuelambda);
      _FA[m]= *K;
  }

}
