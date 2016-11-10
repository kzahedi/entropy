#include "Test.h"

//create values only with random lambda
Test::Test(int colX,int colValY,int rowX,vector<double> lambda,DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY){
  _X= &aX;
  _Y= &aY;
  _systX       = systX;
  _systY       = systY;
  _sizeColValY=colValY;
  __getValX(colX,rowX);
  _sizeColValX=_valX->columns();
  _exact= new IT(colValY,*_valX,*_X,*_Y,systX,systY);
  __setLambdaRand(lambda);
  __getValY(colValY,rowX);
}
//create values only
Test::Test(int colX,int colValY,int rowX,IContainer &indizes, DContainer &lambda,DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY){
  _X= &aX;
  _Y= &aY;
  _systX       = systX;
  _systY       = systY;
  _sizeColValY=colValY;
  __getValX(colX,rowX);
  _sizeColValX=_valX->columns();
  _exact= new IT(colValY,*_valX,*_X,*_Y,systX,systY);
  __setLambda(indizes,lambda);
  __getValY(colValY,rowX);
}
//Angabe der zu vergleichenden Algorithmen, case 0= GIS, 1=SCGIS , 2= GISgp, 3 = SCGISgp
void Test::compareCases(IsParameter param, vector<int>& cases){
  time_t befor;
  time_t after;
  for(int i=0;i<cases.size();i++)
  {
	switch(cases[i]){
	  case 0:
        befor=time(NULL);
        _gisTest= new GIS(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); //1,maxit,konv,test,_timetest,seconds);
	    after=time(NULL);
	    cout <<"GIS: " <<" time: " << difftime(after,befor) <<" distance: " << __KL(0)<< " iterations: " << __getSizeConv(0) << endl;
	    break;
	  case 1:
	    befor=time(NULL);
	    _scgisTest= new SCGIS(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,maxit,konv,test,_timetest,seconds);
	    after=time(NULL);
	    cout <<"SCGIS: " << " time: " << difftime(after,befor) <<" distance: " << __KL(1) <<" iterations: " << __getSizeConv(1) << endl;
        break;
	  case 2:
	    befor=time(NULL);
        _gisgpTest= new GISgp(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,1,0.01,maxit,konv,test,_timetest,seconds);
	    after=time(NULL);
        cout <<"GISgp: " <<" time: " << difftime(after,befor) <<" distance: " << __KL(2) <<" iterations: " << __getSizeConv(2) << endl;
	    break;
	  case 3:
	    befor=time(NULL);
	    _scgisgpTest= new SCGISgp(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,1,0.01,maxit,konv,test,_timetest,seconds);
	    after=time(NULL);
	    cout <<"SCGISgp: " <<" time: " << difftime(after,befor)  <<" distance: " << __KL(3) <<" iterations: " << __getSizeConv(3) << endl;
	    break;
	  default:
	    cout << "default " << endl;
	}
  }
}
Test::~Test()
{
  delete _exact;
  delete _gisTest;
  delete _scgisTest;
  delete _scgisgpTest;
  delete _gisgpTest;
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

double Test::__KL(int i){
  double dist=0;
  double p1=0;
  double q=0;
  double pm1=0;
  assert(abs(i)>=0 && abs(i)<4);
  for(int rowX=0;rowX<pow(_X->rows(),_sizeColValX);rowX++)
  {
    switch(i){
      case 0:
        pm1=_gisTest->propm(rowX);
        break;
      case 1:
        pm1=_scgisTest->propm(rowX);
        break;
      case 2:
        pm1=_gisgpTest->propm(rowX);
        break;
      case 3:
        pm1=_scgisgpTest->propm(rowX);
        break;
      default:
        cout << "default " << endl;
    }
    for(int sizeY=0;sizeY<pow(_Y->rows(),_sizeColValY);sizeY++)
    {
      switch(i){
        case 0:
          p1=_gisTest->propAlphX(rowX,sizeY);
          break;
        case 1:
          p1=_scgisTest->propAlphX(rowX,sizeY);
          break;
        case 2:
          p1=_gisgpTest->propAlphX(rowX,sizeY);
          break;
        case 3:
          p1=_scgisgpTest->propAlphX(rowX,sizeY);
          break;
        default:
          cout << "default " << endl;
      }
      q= _exact->propAlphX(rowX,sizeY);
      if(fabs(p1)<0.00000001){ p1=0.000001;}
      if(fabs( q)<0.00000001){ q =0.000001;}
      dist+=pm1*p1*log(p1/q);
    }
  }
  return dist;
}
//Ausgabe der exakten Wahrscheinlichkeiten
double Test:: getProp(int indexX,int indexY){
	return _exact->propAlphX(indexX,indexY);
}
//fuer die Tests
double Test:: KL(IT *it){
  double dist=0;
  double p1=0;
  double q=0;
  double pm1=0;
  for(int rowX=0;rowX<pow(_X->rows(),_sizeColValX);rowX++){
        pm1=it->propm(rowX);
    for(int sizeY=0;sizeY<pow(_Y->rows(),_sizeColValY);sizeY++)
    {
      p1=it->propAlphX(rowX,sizeY);
      q= _exact->propAlphX(rowX,sizeY);
      if(fabs(p1)<0.00000001){ p1=0.000001;}
      if(fabs( q)<0.00000001){ q =0.000001;}
      dist+=pm1*p1*log(p1/q);
    }
  }
  return dist;
}
//DContainer mit Indizes fuer das lambda und den Wert
void Test:: __setLambda(IContainer &indizes, DContainer &values ){
  assert(indizes.columns() ==3);
  assert(indizes.rows() == values.rows());
  for(int i=0;i< indizes.rows();i++)
  {
    _exact->setFeatureArraylambda(indizes.get(i,0),indizes.get(i,1),indizes.get(i,2),values.get(i,0));
  }
}
//zufaelliges setzen von lambda innerhalb der Werte aus dem Vektor
void Test:: __setLambdaRand(vector<double> lambda){
  srand(time(NULL));
  for(int feat=0;feat<_systX.size();feat++)
  {
    for(int delti=0;delti<pow(_X->rows(),_systX[feat].size());delti++)
    {
      for(int deltj=0; deltj<pow(_Y->rows(),_systY[feat].size());deltj++)
      {
        double r = rand() % lambda.size();
        _exact->setFeatureArraylambda(feat,delti,deltj,lambda[r]);
      }
    }
  }
}

void Test::__getValX(int colX,int rowX)
{
  _valX = new DContainer(rowX,colX);
  for(int i=0;i < rowX; i++)
  {
    for(int j=0;j < colX; j++)
    {
      double z = rand() % _X->rows(); // random indices for the alphabet
      *_valX << (*_X)(z,0);
    }
  }
}

void  Test::__getValY(int colY,int rowX)
{
  srand(time(NULL));
  _valY = new DContainer(rowX,colY);
  for(int i=0;i<rowX;i++)
  {
    double z= (double) rand()/RAND_MAX;
    double s=0;
    int ind=0;
    for(int j=0;j< pow(_Y->rows(),_sizeColValY) && s<z;j++)
    {
      s+= _exact->prop(i,j);;
      ind=j;
    }
    for(int k=0;k<_sizeColValY;k++)
    {
      (*_valY) << _exact->index(ind,false,_sizeColValY)[k] ;
    }
  }
}

DContainer& Test:: getvalX()
{
  return *_valX;
}

DContainer& Test:: getvalY()
{
  return *_valY;
}

int Test::__getSizeConv(int i)
{
  assert(abs(i)>=0 && abs(i)<4);
  switch(i)
  {
    case 0:
      return _gisTest->getsizeconv();
      break;
    case 1:
      return _scgisTest->getsizeconv();
      break;
    case 2:
      return _gisgpTest->getsizeconv();
      break;
    case 3:
      return _scgisgpTest->getsizeconv();
      break;
    default:
      cout << "default " << endl;
  }
  return -1;
}
