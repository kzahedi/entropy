#include "Test.h"

//vergleichswerte, gemessene X,Y und Eingabealphabete
Test::Test(int colX,int colValY, int rowX, vector<double> lambda, DContainer &aX, DContainer &aY, IsParameter param)
{
  _timetest    = param.time;
  _X           = &aX;
  _Y           = &aY;
  _sizeColValY = colValY;
  __getValX(colX,rowX);
  _sizeColValX = _valX->columns();
  _exact       = new IT(colValY,*_valX,*_X,*_Y);
  __setLambdaRand(lambda);
  _alphY       = __getalph(false);
  _alphX       = __getalph(true);
  __getValY(colValY,rowX);
  __comptime(param);
  _case        = 4;
}

Test::Test(int colX,int colValY,int rowX, IContainer &indizes, DContainer &lambda ,DContainer &aX, DContainer &aY, IsParameter param)
{
  _timetest    = param.time;
  _X           = &aX;
  _Y           = &aY;
  _sizeColValY = colValY;
  __getValX(colX,rowX);
  _sizeColValX = _valX->columns();
  _exact       = new IT(colValY,*_valX,*_X,*_Y);
  __setlambda(indizes,lambda);
  _alphY       = __getalph(false);
  _alphX       = __getalph(true);
  __getValY(colValY,rowX);
  __comptime(param);
  _case        = 4;
}

// test if same initial conditions lead to same convergence behaviour
Test::Test(int colX,int colValY,int rowX,DContainer &aX, DContainer &aY, IsParameter param)
{
  _timetest    = param.time;
  _X           = &aX;
  _Y           = &aY;
  _sizeColValY = colValY;
  __getValX(colX,rowX);
  _sizeColValX = _valX->columns();
  _exact= new IT(colValY,*_valX,*_X,*_Y);
  _exact->setFeatureArraylambda(0,0,0,0,1);
  _exact->setFeatureArraylambda(0,0,1,0,4);
  _exact->setFeatureArraylambda(0,0,0,1,5);
  _exact->setFeatureArraylambda(0,0,1,1,0);

  _exact->setFeatureArraylambda(1,0,0,0,0);
  _exact->setFeatureArraylambda(1,0,1,0,2);
  _exact->setFeatureArraylambda(1,0,0,1,1);
  _exact->setFeatureArraylambda(1,0,1,1,3);

  _exact->setFeatureArraylambda(0,1,0,0,3);
  _exact->setFeatureArraylambda(0,1,1,0,2);
  _exact->setFeatureArraylambda(0,1,0,1,0);
  _exact->setFeatureArraylambda(0,1,1,1,3);

  _exact->setFeatureArraylambda(1,1,0,0,2);
  _exact->setFeatureArraylambda(1,1,1,0,1);
  _exact->setFeatureArraylambda(1,1,0,1,0);
  _exact->setFeatureArraylambda(1,1,1,1,3);
  _alphY=__getalph(false);
  _alphX=__getalph(true);
  __getValY(colValY,rowX);
  __comptime(param);
  _case=4;
}

Test::Test(int colX,int colValY,int rowX,vector<double> lambda,DContainer &aX, DContainer &aY){
  _X= &aX;
  _Y= &aY;
  _sizeColValY=colValY;
  __getValX(colX,rowX);
  _sizeColValX=_valX->columns();
  _exact= new IT(colValY,*_valX,*_X,*_Y);
  __setLambdaRand(lambda);
  _alphY=__getalph(false);
  _alphX=__getalph(true);
  __getValY(colValY,rowX);
  _case=5;
}

Test::Test(int colX, int colValY, int rowX, vector<double> lambda, DContainer &aX, DContainer &aY, IsParameter param, int i) // int maxit, double konv, bool time,bool test,int seconds,int i){
{
  assert(abs(i)>=0 && abs(i)<4);
  _case        = i;
  _timetest    = param.time;
  _X           = &aX;
  _Y           = &aY;
  _sizeColValY = colValY;
  __getValX(colX,rowX);
  _sizeColValX = _valX->columns();
  _exact       = new IT(colValY,*_valX,*_X,*_Y);
  __setLambdaRand(lambda);
  _alphY       = __getalph(false);
  _alphX       = __getalph(true);
  __getValY(colValY,rowX);
  switch(i){
    case 0:
      _gisTest=new GIS(*_valX, *_valY, *_X, *_Y, param); //1,maxit,konv,test,_timetest,seconds);
      break;
    case 1:
      _scgisTest=new SCGIS(*_valX,*_valY,*_X,*_Y, param); // 1,maxit,konv,test,_timetest,seconds);
      break;
    case 2:
      _gisgpTest=new GISgp(*_valX,*_valY,*_X,*_Y, param); // 1,1,0.01,maxit,konv,test,_timetest,seconds);
      break;
    case 3:
      _scgisgpTest= new SCGISgp(*_valX,*_valY,*_X,*_Y,param); // 1,1,0.01,maxit,konv,test,_timetest,seconds);
      break;
    default:
      cout << "default " << endl;
  }
}
Test::~Test()
{
  _timediff.clear();
  for(int i=0;i<_alphX.size();i++){
    _alphX[i].clear();
  }
  _alphX.clear();

  for(int i=0;i<_alphY.size();i++){
    _alphY[i].clear();
  }
  _alphY.clear();
  //expizit aufrufen?
  delete _exact;
  delete _gisTest;
  delete _scgisTest;
  delete _scgisgpTest;
  delete _gisgpTest;
}
//GIS und SCGIS ausfuehren mit Zeitmessung
void Test::__comptime(IsParameter param)
{
  time_t befor1;
  time_t befor2;
  time_t after1;
  time_t after2;
  time_t befor3;
  time_t after3;
  time_t befor4;
  time_t after4;
  befor1=time(NULL);
  _gisTest=new GIS(*_valX,*_valY,*_X,*_Y, param); // 1,maxit,konv,false,_timetest,seconds);
  after1=time(NULL);
  _timediff.push_back(difftime(after1,befor1));
  befor2=time(NULL);
  _scgisTest=new SCGIS(*_valX,*_valY,*_X,*_Y, param); // 1,maxit,konv,false,_timetest,seconds);
  after2=time(NULL);
  _timediff.push_back(difftime(after2,befor2));
  befor3=time(NULL);
  _gisgpTest=new GISgp(*_valX,*_valY,*_X,*_Y, param); // 1,1,0.01,maxit,konv,false,_timetest,seconds);
  after3=time(NULL);
  _timediff.push_back(difftime(after3,befor3));
  befor4=time(NULL);
  _scgisgpTest= new SCGISgp(*_valX,*_valY,*_X,*_Y, param); // 1,1,0.01,maxit,konv,false,_timetest,seconds);
  after4=time(NULL);
  _timediff.push_back(difftime(after4,befor4));
}
//Ausgabe
void Test:: comparison()
{
  cout << "Comparison: " << endl;
  cout << endl;
  cout << "time:  GIS: " << _timediff[0] << "s SCGIS: " << _timediff[1] << "s GIS smoothed: " << _timediff[2]<<  "s SCGIS smoothed: "  << _timediff[3] << "s" << endl;
  cout << endl;
  vector<double> kl = KL();
  cout << "KL-distance: GIS: " << kl[0] <<  " SCGIS: " << kl[1] <<" GIS smoothed: "<< kl[2] << " SCGIS smoothed: " << kl[3] <<endl;
  cout << endl;
  if(_timetest){
    cout<< "Iterations: GIS: " << _gisTest->getIterations() << " SCGIS: " << _scgisTest->getIterations()<< " GIS smoothed: " << _gisgpTest->getIterations() << " SCGIS smoothed: " << _scgisgpTest->getIterations() <<   endl;
  }
  // cout << "lambda: " << endl;
  // cout << "vergleichswerte" << endl;
  // cout << _exact->getFeatureArraylambda(0,0,0,0) <<endl;
  // cout << _exact->getFeatureArraylambda(0,0,1,0) <<endl;
  // cout << _exact->getFeatureArraylambda(0,0,0,1) <<endl;
  // cout << _exact->getFeatureArraylambda(0,0,1,1) <<endl;
  // cout << endl;

  // cout << _exact->getFeatureArraylambda(1,0,0,0) <<endl;
  // cout << _exact->getFeatureArraylambda(1,0,1,0) <<endl;
  // cout << _exact->getFeatureArraylambda(1,0,0,1) <<endl;
  // cout << _exact->getFeatureArraylambda(1,0,1,1) <<endl;
  // cout << endl;
  // cout << _exact->getFeatureArraylambda(0,1,0,0) <<endl;
  // cout << _exact->getFeatureArraylambda(0,1,1,0) <<endl;
  // cout << _exact->getFeatureArraylambda(0,1,0,1) <<endl;
  // cout << _exact->getFeatureArraylambda(0,1,1,1) <<endl;
  // cout << endl;
  // cout << _exact->getFeatureArraylambda(1,1,0,0) <<endl;
  // cout << _exact->getFeatureArraylambda(1,1,1,0) <<endl;
  // cout << _exact->getFeatureArraylambda(1,1,0,1) <<endl;
  // cout << _exact->getFeatureArraylambda(1,1,1,1) <<endl;

  // cout << "GIS" << endl;
  // cout << _gisTest->getFeatureArraylambda(0,0,0,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(0,0,1,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(0,0,0,1) <<endl;
  // cout << _gisTest->getFeatureArraylambda(0,0,1,1) <<endl;
  // cout << endl;

  // cout << _gisTest->getFeatureArraylambda(1,0,0,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(1,0,1,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(1,0,0,1) <<endl;
  // cout << _gisTest->getFeatureArraylambda(1,0,1,1) <<endl;
  // cout << endl;
  // cout << _gisTest->getFeatureArraylambda(0,1,0,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(0,1,1,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(0,1,0,1) <<endl;
  // cout << _gisTest->getFeatureArraylambda(0,1,1,1) <<endl;
  // cout << endl;
  // cout << _gisTest->getFeatureArraylambda(1,1,0,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(1,1,1,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(1,1,0,1) <<endl;
  // cout << _gisTest->getFeatureArraylambda(1,1,1,1) <<endl;
  // cout << "GIS smoothed " << endl;
  // cout << _gisgpTest->getFeatureArraylambda(0,0,0,0) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(0,0,1,0) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(0,0,0,1) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(0,0,1,1) <<endl;
  // cout << endl;
  // cout << _gisgpTest->getFeatureArraylambda(1,0,0,0) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(1,0,1,0) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(1,0,0,1) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(1,0,1,1) <<endl;
  // cout << endl;
  // cout << _gisgpTest->getFeatureArraylambda(0,1,0,0) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(0,1,1,0) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(0,1,0,1) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(0,1,1,1) <<endl;
  // cout << endl;
  // cout << _gisgpTest->getFeatureArraylambda(1,1,0,0) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(1,1,1,0) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(1,1,0,1) <<endl;
  // cout << _gisgpTest->getFeatureArraylambda(1,1,1,1) <<endl;

  // cout << "alphY :" << endl;
  // for(int i=0;i<_alphY.size();i++){
    // for(int j=0;j<_alphY[i].size();j++){
      // cout << _alphY[i][j] << " " ;
    // }
    // cout << endl;
  // }
  // cout << "alphX :" << endl;
  // for(int i=0;i<_alphX.size();i++){
    // for(int j=0;j<_alphX[i].size();j++){
      // cout << _alphX[i][j] << " " ;
    // }
    // cout << endl;
  // }
  // cout << _alphX.size() << endl;

}
//Abstand
vector<double> Test:: KL(){
  vector<double> dist(4);
  double p1=0;
  double p2=0;
  double p3=0;
  double p4=0;
  double q=0;
  double pm1=0;
  double pm2=0;
  double pm3=0;
  double pm4=0;
  for(int rowX=0;rowX<  _alphX.size() ;rowX++)
  {
    pm1=_gisTest->propm(_alphX,rowX,_alphY);
    pm2=_scgisTest->propm(_alphX,rowX,_alphY);
    pm3=_gisgpTest->propm(_alphX,rowX,_alphY);
    pm4=_scgisgpTest->propm(_alphX,rowX,_alphY);
    for(int sizeY=0;sizeY<_alphY.size();sizeY++)
    {
      p1=_gisTest->prop(rowX,_alphY,sizeY);
      p2=_scgisTest->prop(rowX,_alphY,sizeY);
      p3=_gisgpTest->prop(rowX,_alphY,sizeY);
      p4=_scgisgpTest->prop(rowX,_alphY,sizeY);
      q= _exact->prop(rowX,_alphY,sizeY);
      if(fabs(p1)<0.00000001){ p1=0.000001;}
      if(fabs(p2)<0.00000001){ p2=0.000001;}
      if(fabs(p3)<0.00000001){ p3=0.000001;}
      if(fabs( q)<0.00000001){ q =0.000001;}
      dist[0]+=pm1*p1*log(p1/q);
      dist[1]+=pm2*p2*log(p2/q);
      dist[2]+=pm3*p3*log(p3/q);
      dist[3]+=pm4*p4*log(p4/q);
    }
  }
  return dist;
}
double Test::KL1(){
  double dist=0;
  double p1=0;
  double q=0;
  double pm1=0;
  assert(abs(_case)>=0 && abs(_case)<4);
  for(int rowX=0;rowX<_alphX.size();rowX++){
    switch(_case){
      case 0:
        pm1=_gisTest->propm(_alphX,rowX,_alphY);
        break;
      case 1:
        pm1=_scgisTest->propm(_alphX,rowX,_alphY);
        break;
      case 2:
        pm1=_gisgpTest->propm(_alphX,rowX,_alphY);
        break;
      case 3:
        pm1=_scgisgpTest->propm(_alphX,rowX,_alphY);
        break;
      default:
        cout << "default " << endl;
    }
    for(int sizeY=0;sizeY<_alphY.size();sizeY++){
      switch(_case){
        case 0:
          p1=_gisTest->prop(rowX,_alphY,sizeY);
          break;
        case 1:
          p1=_scgisTest->prop(rowX,_alphY,sizeY);
          break;
        case 2:
          p1=_gisgpTest->prop(rowX,_alphY,sizeY);
          break;
        case 3:
          p1=_scgisgpTest->prop(rowX,_alphY,sizeY);
          break;
        default:
          cout << "default " << endl;
      }
      q= _exact->prop(rowX,_alphY,sizeY);
      if(fabs(p1)<0.00000001){ p1=0.000001;}
      if(fabs( q)<0.00000001){ q =0.000001;}
      dist+=pm1*p1*log(p1/q);
    }
  }
  return dist;
}
// true fuer x, false fuer y
vector<vector<double> > Test::__getalph(bool valX)
{
  int rowsAlph;
  int colval;
  if(valX)
  {
    rowsAlph = _X->rows();
    colval   = _sizeColValX;
  }
  else
  {
    rowsAlph = _Y->rows();
    colval   = _sizeColValY;
  }
  vector<double> fill(0);
  vector<vector<double> > Z(pow(rowsAlph,colval),vector<double>(colval));
  _index=0;
  __fill(fill,0,Z,valX,rowsAlph,colval);
  return Z;
}

void Test::__fill(vector<double> fill, int i, vector<vector<double> > &Z,bool valX,int rowsAlph,int colval)
{
  if(colval==i)
  {
    for(int l=0;l<colval;l++)
    {
      Z[_index][l]=fill[l];
    }
    _index++;
    fill.pop_back();
    i--;
  }
  else
  {
    for(int x=0;x<rowsAlph;x++)
    {
      if(valX){
        fill.push_back(_X->get(x,0));
      }
      else{
        fill.push_back(_Y->get(x,0));
      }
      __fill(fill,i+1,Z,valX,rowsAlph,colval);
      fill.pop_back();
    }
  }
}
//DContainer mit Indizes fuer das lambda und den Wert
void Test:: __setlambda(IContainer &indizes, DContainer &values ){
  assert(indizes.columns() ==4);
  assert(indizes.rows() == values.rows());
  for(int i=0;i< indizes.rows();i++){
    _exact->setFeatureArraylambda(indizes.get(i,0),indizes.get(i,1),indizes.get(i,2),indizes.get(i,3),values.get(i,0));
  }
}
void Test:: __setLambdaRand(vector<double> lambda){
  srand(time(NULL));
  for(int feati=0;feati<_sizeColValX;feati++){
    for(int featj=0; featj<_sizeColValY;featj++){
      for(int delti=0;delti<_X->rows();delti++){
        for(int deltj=0; deltj<_Y->rows();deltj++){
          double r = rand() % lambda.size();
          _exact->setFeatureArraylambda(feati,featj,delti,deltj,lambda[r]);
        }
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
  double** prop;
  prop=new double*[rowX];
  for(int i=0;i<rowX;i++){
    prop[i]=new double[_alphY.size()];
    for(int j=0;j<_alphY.size();j++){
      prop[i][j]=0;
    }
  }
  for(int i=0;i<rowX;i++){
    for(int propi=0;propi<_alphY.size();propi++){
      double l=_exact->prop(i,_alphY,propi);
      prop[i][propi]=_exact->prop(i,_alphY,propi);
    }
  }
  _valY = new DContainer(rowX,colY);
  for(int i=0;i<rowX;i++){
    double z= (double) rand()/RAND_MAX;
    double s=0;
    int ind=0;
    for(int j=0;j<_alphY.size() && s<z;j++){
      s+= prop[i][j];
      ind=j;
    }
    for(int k=0;k<_sizeColValY;k++){
      (*_valY) << _alphY[ind][k];
    }
  }
  for(int i=0;i<rowX;i++){
    delete [] prop[i];
  }
  delete [] prop;
}

DContainer& Test:: getvalX()
{
  return *_valX;
}

DContainer& Test:: getvalY()
{
  return *_valY;
}

double Test::prop(int feati,int featj,double valX,double valY)
{
  assert(abs(_case)>=0 && abs(_case)<4);
  switch(_case)
  {
    case 0:
      return _gisTest->prop(feati,featj,valX,valY);
      break;
    case 1:
      return _scgisTest->prop(feati,featj,valX,valY);
      break;
    case 2:
      return _gisgpTest->prop(feati,featj,valX,valY);
      break;
    case 3:
      return _scgisgpTest->prop(feati,featj,valX,valY);
      break;
    default:
      cout << "default " << endl;
  }
  return 0.0;
}

double Test:: getconv(int ind)
{
  assert(abs(_case)>=0 && abs(_case)<4);
  switch(_case){
    case 0:
      return _gisTest->getconv(ind);
      break;
    case 1:
      return _scgisTest->getconv(ind);
      break;
    case 2:
      return _gisgpTest->getconv(ind);
      break;
    case 3:
      return _scgisTest->getconv(ind);
      break;
    default:
      cout << "default " << endl;
  }
  return 0.0;
}

int Test::getsizeconv()
{
  assert(abs(_case)>=0 && abs(_case)<4);
  switch(_case)
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
      return _scgisTest->getsizeconv();
      break;
    default:
      cout << "default " << endl;
  }
  return -1;
}

