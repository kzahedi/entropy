#include "Test.h"
/*
//vergleichswerte, gemessene X,Y und Eingabealphabete
Test::Test(int colX,int colValY, int rowX, vector<double> lambda, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
{
  _timetest    = param.time;
  _X           = &aX;
  _Y           = &aY;
  _systX       = systX;
  _systY       = systY;
  _sizeColValY = colValY;
  __getValX(colX,rowX);
  _sizeColValX = _valX->columns();
  _exact       = new IT(colValY,*_valX,*_X,*_Y,_systX,_systY);
  //__setLambdaRand(lambda);

	 _exact->setFeatureArraylambda(0,0,0,1);
	 _exact->setFeatureArraylambda(0,1,0,4);
	 _exact->setFeatureArraylambda(0,0,1,5);
	 _exact->setFeatureArraylambda(0,1,1,0);

	 _exact->setFeatureArraylambda(1,0,0,0);
	 _exact->setFeatureArraylambda(1,1,0,2);
	 _exact->setFeatureArraylambda(1,0,1,1);
	 _exact->setFeatureArraylambda(1,1,1,3);

	 _exact->setFeatureArraylambda(2,0,0,3);
	 _exact->setFeatureArraylambda(2,1,0,2);
	 _exact->setFeatureArraylambda(2,0,1,0);
	 _exact->setFeatureArraylambda(2,1,1,3);

	 _exact->setFeatureArraylambda(3,0,0,2);
	 _exact->setFeatureArraylambda(3,1,0,1);
	 _exact->setFeatureArraylambda(3,0,1,0);
	 _exact->setFeatureArraylambda(3,1,1,3);
  __getValY(colValY,rowX);
  __comptime(param);
  _case        = 4;
}

Test::Test(int colX,int colValY,int rowX, IContainer &indizes, DContainer &lambda ,DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
{
  _timetest    = param.time;
  _X           = &aX;
  _Y           = &aY;
  _systX       = systX;
  _systY       = systY;
  _sizeColValY = colValY;
  __getValX(colX,rowX);
  _sizeColValX = _valX->columns();
  _exact       = new IT(colValY,*_valX,*_X,*_Y,systX,systY);
  __setlambda(indizes,lambda);
  __getValY(colValY,rowX);
  __comptime(param);
  _case        = 4;
}
*/
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
  _case=5;
}
void Test::compareCases( IsParameter param, vector<int>& cases){
  time_t befor;
  time_t after;
  for(int i=0;i<cases.size();i++){
	  switch(cases[i]){
	      case 0:
	        befor=time(NULL);
	        _gisTest= new GIS(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); //1,maxit,konv,test,_timetest,seconds);
	        after=time(NULL);
	        cout <<"GIS: ";
	        cout <<" time: " << difftime(after,befor);
	        cout <<" distance: " << __KL(0);
	        cout <<" iterations: " << __getSizeConv(0) << endl;
	        break;
	      case 1:
	        befor=time(NULL);
	        _scgisTest= new SCGIS(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,maxit,konv,test,_timetest,seconds);
	        after=time(NULL);
	        cout <<"SCGIS: ";
	        cout <<" time: " << difftime(after,befor);
	        cout <<" distance: " << __KL(1);
	        cout <<" iterations: " << __getSizeConv(1) << endl;
	        break;
	      case 2:
	    	  cout << " gisgp Test" << endl;
	        befor=time(NULL);
	        _gisgpTest= new GISgp(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,1,0.01,maxit,konv,test,_timetest,seconds);
	        after=time(NULL);
	        cout <<"GISgp: ";
	        cout <<" time: " << difftime(after,befor);
	        cout <<" distance: " << __KL(2);
	        cout <<" iterations: " << __getSizeConv(2) << endl;
	        break;
	      case 3:
	    	  cout << " vor scgisgptest" << endl;
	      	befor=time(NULL);
	        _scgisgpTest= new SCGISgp(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,1,0.01,maxit,konv,test,_timetest,seconds);
	        after=time(NULL);
	        cout <<"SCGISgp: ";
	        cout <<" time: " << difftime(after,befor);
	        cout <<" distance: " << __KL(3);
	        cout <<" iterations: " << __getSizeConv(3) << endl;
	        break;
	      default:
	        cout << "default " << endl;
	    }
  }
}
/*
Test::Test(int colX, int colValY, int rowX, vector<double> lambda, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param, int i)
{
  assert(abs(i)>=0 && abs(i)<4);
  _case        = i;
  _timetest    = param.time;
  _X           = &aX;
  _Y           = &aY;
  _systX       = systX;
  _systY       = systY;
  _sizeColValY = colValY;
  __getValX(colX,rowX);
  _sizeColValX = _valX->columns();
  _exact       = new IT(colValY,*_valX,*_X,*_Y,_systX,_systY);
  //__setLambdaRand(lambda);
	 _exact->setFeatureArraylambda(0,0,0,1);
	 _exact->setFeatureArraylambda(0,1,0,4);
	 _exact->setFeatureArraylambda(0,0,1,5);
	 _exact->setFeatureArraylambda(0,1,1,0);

	 _exact->setFeatureArraylambda(1,0,0,0);
	 _exact->setFeatureArraylambda(1,1,0,2);
	 _exact->setFeatureArraylambda(1,0,1,1);
	 _exact->setFeatureArraylambda(1,1,1,3);

	 _exact->setFeatureArraylambda(2,0,0,3);
	 _exact->setFeatureArraylambda(2,1,0,2);
	 _exact->setFeatureArraylambda(2,0,1,0);
	 _exact->setFeatureArraylambda(2,1,1,3);

	 _exact->setFeatureArraylambda(3,0,0,2);
	 _exact->setFeatureArraylambda(3,1,0,1);
	 _exact->setFeatureArraylambda(3,0,1,0);
	 _exact->setFeatureArraylambda(3,1,1,3);
	  time_t befor;
	  time_t after;
  __getValY(colValY,rowX);
  switch(i){
    case 0:
      befor=time(NULL);
      _gisTest= new GIS(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); //1,maxit,konv,test,_timetest,seconds);
      after=time(NULL);
      cout << difftime(after,befor);
      break;
    case 1:
      befor=time(NULL);
      _scgisTest= new SCGIS(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,maxit,konv,test,_timetest,seconds);
      after=time(NULL);
      cout << difftime(after,befor);
      break;
    case 2:
      befor=time(NULL);
      _gisgpTest= new GISgp(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,1,0.01,maxit,konv,test,_timetest,seconds);
      after=time(NULL);
      cout << difftime(after,befor);
      break;
    case 3:
    	befor=time(NULL);
      _scgisgpTest= new SCGISgp(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,1,0.01,maxit,konv,test,_timetest,seconds);
      after=time(NULL);
      cout<< difftime(after,befor);
      break;
    default:
      cout << "default " << endl;
  }
} */
Test::~Test()
{
  _timediff.clear();
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
/*
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
  cout << " vor GIS " << endl;
  _gisTest=new GIS(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,maxit,konv,false,_timetest,seconds);
  after1=time(NULL);
  _timediff.push_back(difftime(after1,befor1));
  befor2=time(NULL);
  cout << " vor SCGIS " << endl;
  _scgisTest=new SCGIS(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,maxit,konv,false,_timetest,seconds);
  after2=time(NULL);
  _timediff.push_back(difftime(after2,befor2));
  befor3=time(NULL);
  cout << " vor GISgp " << endl;
  _gisgpTest=new GISgp(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,1,0.01,maxit,konv,false,_timetest,seconds);
  after3=time(NULL);
  _timediff.push_back(difftime(after3,befor3));
  befor4=time(NULL);
  cout << " vor SCGISgp " << endl;
  _scgisgpTest= new SCGISgp(*_valX,*_valY,*_X,*_Y,_systX,_systY, param); // 1,1,0.01,maxit,konv,false,_timetest,seconds);
  after4=time(NULL);
  _timediff.push_back(difftime(after4,befor4));
}
//Ausgabe
void Test:: comparison()
{
  cout << "Comparison: " << endl;
  cout << endl;
  cout << "time:  GIS: " << _timediff[0]  << "s SCGIS: " << _timediff[1] << "s GIS smoothed: " << _timediff[2] <<  "s SCGIS smoothed: "  << _timediff[3] << "s" << endl;
  cout << endl;
  vector<double> kl = KL();
  cout << "KL-distance: GIS: " << kl[0] <<  " SCGIS: " << kl[1] <<" GIS smoothed: "<< kl[2] << " SCGIS smoothed: " << kl[3] <<endl;
  cout << endl;
  cout << KL1(0) << " " << KL1(1) <<" " <<  KL1(2) << " " << KL1(3) << endl;
  //if(_timetest){
    cout<< "Iterations: GIS: " << _gisTest->getIterations() << " SCGIS: " << _scgisTest->getIterations() << " GIS smoothed: " << _gisgpTest->getIterations() << " SCGIS smoothed: " << _scgisgpTest->getIterations() <<   endl;
//  }

   cout << "lambda: " << endl;
   cout << "vergleichswerte" << endl;
   cout << _exact->getFeatureArraylambda(0,0,0) <<endl;
   cout << _exact->getFeatureArraylambda(0,1,0) <<endl;
  // cout << _exact->getFeatureArraylambda(0,2,0) <<endl;
  // cout << _exact->getFeatureArraylambda(0,3,0) <<endl;
   cout << _exact->getFeatureArraylambda(0,0,1) <<endl;
   cout << _exact->getFeatureArraylambda(0,1,1) <<endl;
 //  cout << _exact->getFeatureArraylambda(0,2,1) <<endl;
 //  cout << _exact->getFeatureArraylambda(0,3,1) <<endl;
   cout << endl;

   cout << _exact->getFeatureArraylambda(1,0,0) <<endl;
   cout << _exact->getFeatureArraylambda(1,1,0) <<endl;
   cout << _exact->getFeatureArraylambda(1,0,1) <<endl;
   cout << _exact->getFeatureArraylambda(1,1,1) <<endl;
   cout << endl;
   cout << _exact->getFeatureArraylambda(2,0,0) <<endl;
   cout << _exact->getFeatureArraylambda(2,1,0) <<endl;
   cout << _exact->getFeatureArraylambda(2,0,1) <<endl;
   cout << _exact->getFeatureArraylambda(2,1,1) <<endl;
   cout << endl;
   cout << _exact->getFeatureArraylambda(3,0,0) <<endl;
   cout << _exact->getFeatureArraylambda(3,1,0) <<endl;
   cout << _exact->getFeatureArraylambda(3,0,1) <<endl;
   cout << _exact->getFeatureArraylambda(3,1,1) <<endl;

   cout << "GIS" << endl;
   cout << _gisTest->getFeatureArraylambda(0,0,0) <<endl;
   cout << _gisTest->getFeatureArraylambda(0,1,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(0,2,0) <<endl;
  // cout << _gisTest->getFeatureArraylambda(0,3,0) <<endl;
   cout << _gisTest->getFeatureArraylambda(0,0,1) <<endl;
   cout << _gisTest->getFeatureArraylambda(0,1,1) <<endl;
 //  cout << _gisTest->getFeatureArraylambda(0,2,1) <<endl;
 //  cout << _gisTest->getFeatureArraylambda(0,3,1) <<endl;
   cout << endl;

   cout << _gisTest->getFeatureArraylambda(1,0,0) <<endl;
   cout << _gisTest->getFeatureArraylambda(1,1,0) <<endl;
   cout << _gisTest->getFeatureArraylambda(1,0,1) <<endl;
   cout << _gisTest->getFeatureArraylambda(1,1,1) <<endl;
   cout << endl;

   cout << _gisTest->getFeatureArraylambda(2,0,0) <<endl;
   cout << _gisTest->getFeatureArraylambda(2,1,0) <<endl;
   cout << _gisTest->getFeatureArraylambda(2,0,1) <<endl;
   cout << _gisTest->getFeatureArraylambda(2,1,1) <<endl;
   cout << endl;
   cout << endl;
   cout << _gisTest->getFeatureArraylambda(3,0,0) <<endl;
   cout << _gisTest->getFeatureArraylambda(3,1,0) <<endl;
   cout << _gisTest->getFeatureArraylambda(3,0,1) <<endl;
   cout << _gisTest->getFeatureArraylambda(3,1,1) <<endl;

   cout << "SCGIS" << endl;
   cout << _scgisTest->getFeatureArraylambda(0,0,0) <<endl;
   cout << _scgisTest->getFeatureArraylambda(0,1,0) <<endl;
 //  cout << _scgisTest->getFeatureArraylambda(0,2,0) <<endl;
 //  cout << _scgisTest->getFeatureArraylambda(0,3,0) <<endl;
   cout << _scgisTest->getFeatureArraylambda(0,0,1) <<endl;
   cout << _scgisTest->getFeatureArraylambda(0,1,1) <<endl;
 //  cout << _scgisTest->getFeatureArraylambda(0,2,1) <<endl;
 //  cout << _scgisTest->getFeatureArraylambda(0,3,1) <<endl;
   cout << endl;

   cout << _scgisTest->getFeatureArraylambda(1,0,0) <<endl;
   cout << _scgisTest->getFeatureArraylambda(1,1,0) <<endl;
   cout << _scgisTest->getFeatureArraylambda(1,0,1) <<endl;
   cout << _scgisTest->getFeatureArraylambda(1,1,1) <<endl;
   cout << endl;

   cout << _scgisTest->getFeatureArraylambda(2,0,0) <<endl;
   cout << _scgisTest->getFeatureArraylambda(2,1,0) <<endl;
   cout << _scgisTest->getFeatureArraylambda(2,0,1) <<endl;
   cout << _scgisTest->getFeatureArraylambda(2,1,1) <<endl;
   cout << endl;

   cout << _scgisTest->getFeatureArraylambda(3,0,0) <<endl;
   cout << _scgisTest->getFeatureArraylambda(3,1,0) <<endl;
   cout << _scgisTest->getFeatureArraylambda(3,0,1) <<endl;
   cout << _scgisTest->getFeatureArraylambda(3,1,1) <<endl;
   cout << endl;
/*
   cout << "GIS smoothed " << endl;
   cout << _gisgpTest->getFeatureArraylambda(0,0,0) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(0,1,0) <<endl;
 //  cout << _gisgpTest->getFeatureArraylambda(0,2,0) <<endl;
 //  cout << _gisgpTest->getFeatureArraylambda(0,3,0) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(0,0,1) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(0,1,1) <<endl;
//   cout << _gisgpTest->getFeatureArraylambda(0,2,1) <<endl;
 //  cout << _gisgpTest->getFeatureArraylambda(0,3,1) <<endl;
   cout << endl;

   cout << _gisgpTest->getFeatureArraylambda(1,0,0) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(1,1,0) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(1,0,1) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(1,1,1) <<endl;
   cout << endl;

   cout << _gisgpTest->getFeatureArraylambda(2,0,0) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(2,1,0) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(2,0,1) <<endl;
   cout << _gisgpTest->getFeatureArraylambda(2,1,1) <<endl;
   cout << endl;

 //  cout << _gisgpTest->getFeatureArraylambda(3,0,0) <<endl;
 //  cout << _gisgpTest->getFeatureArraylambda(1,1,1,0) <<endl;
 //  cout << _gisgpTest->getFeatureArraylambda(1,1,0,1) <<endl;
 //  cout << _gisgpTest->getFeatureArraylambda(1,1,1,1) <<endl;


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
  for(int rowX=0;rowX< pow(_X->rows(),_sizeColValX) ;rowX++)
  {
    pm1=_gisTest->propm(rowX);
    pm2=_scgisTest->propm(rowX);
    pm3=_gisgpTest->propm(rowX);
    pm4=_scgisgpTest->propm(rowX);
    for(int sizeY=0;sizeY<pow(_Y->rows(),_sizeColValY);sizeY++)
    {
      p1=_gisTest->prop(rowX,sizeY);
      p2=_scgisTest->prop(rowX,sizeY);
      p3=_gisgpTest->prop(rowX,sizeY);
      p4=_scgisgpTest->prop(rowX,sizeY);
      q= _exact->prop(rowX,sizeY);
      if(fabs(p1)<0.00000001){ p1=0.000001;}
      if(fabs(p2)<0.00000001){ p2=0.000001;}
      if(fabs(p3)<0.00000001){ p3=0.000001;}
      if(fabs(p4)<0.00000001){ p4=0.000001;}
      if(fabs( q)<0.00000001){ q =0.000001;}
      dist[0]+=pm1*p1*log(p1/q);
      dist[1]+=pm2*p2*log(p2/q);
      dist[2]+=pm3*p3*log(p3/q);
      dist[3]+=pm4*p4*log(p4/q);
    }
  }
  //cout << dist[0] << " " << dist[1] <<" " << dist[2] << " " << dist[3] << endl;
  return dist;
} */
double Test::__KL(int i){
  double dist=0;
  double p1=0;
  double q=0;
  double pm1=0;
  assert(abs(i)>=0 && abs(i)<4);
  for(int rowX=0;rowX<pow(_X->rows(),_sizeColValX);rowX++){
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
          p1=_gisTest->prop(rowX,sizeY);
          break;
        case 1:
          p1=_scgisTest->prop(rowX,sizeY);
          break;
        case 2:
          p1=_gisgpTest->prop(rowX,sizeY);
          break;
        case 3:
          p1=_scgisgpTest->prop(rowX,sizeY);
          break;
        default:
          cout << "default " << endl;
      }
      q= _exact->prop(rowX,sizeY);
      if(fabs(p1)<0.00000001){ p1=0.000001;}
      if(fabs( q)<0.00000001){ q =0.000001;}
      dist+=pm1*p1*log(p1/q);
    }
  }
  return dist;
}

//DContainer mit Indizes fuer das lambda und den Wert
void Test:: __setlambda(IContainer &indizes, DContainer &values ){
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
  for(int i=0;i<rowX;i++){
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
/*
double Test::prop(int indexX, int indexY)
{
  assert(abs(_case)>=0 && abs(_case)<4);
  switch(_case)
  {
    case 0:
      return _gisTest->prop( indexX, indexY);
      break;
    case 1:
      return _scgisTest->prop( indexX, indexY);
      break;
    case 2:
      return _gisgpTest->prop(indexX, indexY);
      break;
    case 3:
      return _scgisgpTest->prop(indexX, indexY);
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
      return _scgisgpTest->getconv(ind);
      break;
    default:
      cout << "default " << endl;
  }
  return 0.0;
}
*/
int Test::__getSizeConv(int i)
{
  assert(abs(i)>=0 && abs(i)<4);
  switch(i)
  {
    case 0:
     /* cout << "GIS" << endl;
      cout << _gisTest->getFeatureArraylambda(0,0,0) <<endl;
      cout << _gisTest->getFeatureArraylambda(0,1,0) <<endl;
     // cout << _gisTest->getFeatureArraylambda(0,2,0) <<endl;
     // cout << _gisTest->getFeatureArraylambda(0,3,0) <<endl;
      cout << _gisTest->getFeatureArraylambda(0,0,1) <<endl;
      cout << _gisTest->getFeatureArraylambda(0,1,1) <<endl;
    //  cout << _gisTest->getFeatureArraylambda(0,2,1) <<endl;
    //  cout << _gisTest->getFeatureArraylambda(0,3,1) <<endl;
      cout << endl;

      cout << _gisTest->getFeatureArraylambda(1,0,0) <<endl;
      cout << _gisTest->getFeatureArraylambda(1,1,0) <<endl;
      cout << _gisTest->getFeatureArraylambda(1,0,1) <<endl;
      cout << _gisTest->getFeatureArraylambda(1,1,1) <<endl;
      cout << endl;

      cout << _gisTest->getFeatureArraylambda(2,0,0) <<endl;
      cout << _gisTest->getFeatureArraylambda(2,1,0) <<endl;
      cout << _gisTest->getFeatureArraylambda(2,0,1) <<endl;
      cout << _gisTest->getFeatureArraylambda(2,1,1) <<endl;
      cout << endl;
      cout << endl;
      cout << _gisTest->getFeatureArraylambda(3,0,0) <<endl;
      cout << _gisTest->getFeatureArraylambda(3,1,0) <<endl;
      cout << _gisTest->getFeatureArraylambda(3,0,1) <<endl;
      cout << _gisTest->getFeatureArraylambda(3,1,1) <<endl; */
      return _gisTest->getsizeconv();
      break;
    case 1:

     /* cout << "SCGIS" << endl;
      cout << _scgisTest->getFeatureArraylambda(0,0,0) <<endl;
      cout << _scgisTest->getFeatureArraylambda(0,1,0) <<endl;
    //  cout << _scgisTest->getFeatureArraylambda(0,2,0) <<endl;
    //  cout << _scgisTest->getFeatureArraylambda(0,3,0) <<endl;
      cout << _scgisTest->getFeatureArraylambda(0,0,1) <<endl;
      cout << _scgisTest->getFeatureArraylambda(0,1,1) <<endl;
    //  cout << _scgisTest->getFeatureArraylambda(0,2,1) <<endl;
    //  cout << _scgisTest->getFeatureArraylambda(0,3,1) <<endl;
      cout << endl;

      cout << _scgisTest->getFeatureArraylambda(1,0,0) <<endl;
      cout << _scgisTest->getFeatureArraylambda(1,1,0) <<endl;
      cout << _scgisTest->getFeatureArraylambda(1,0,1) <<endl;
      cout << _scgisTest->getFeatureArraylambda(1,1,1) <<endl;
      cout << endl;

      cout << _scgisTest->getFeatureArraylambda(2,0,0) <<endl;
      cout << _scgisTest->getFeatureArraylambda(2,1,0) <<endl;
      cout << _scgisTest->getFeatureArraylambda(2,0,1) <<endl;
      cout << _scgisTest->getFeatureArraylambda(2,1,1) <<endl;
      cout << endl;

      cout << _scgisTest->getFeatureArraylambda(3,0,0) <<endl;
      cout << _scgisTest->getFeatureArraylambda(3,1,0) <<endl;
      cout << _scgisTest->getFeatureArraylambda(3,0,1) <<endl;
      cout << _scgisTest->getFeatureArraylambda(3,1,1) <<endl;
      cout << endl; */
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
/*
//einzelausgaben fuer it_test
vector<double> Test::propAll(int indexX, int indexY){
	assert(_case==4);
	vector<double> prop(4);
    prop[0]=_gisTest->prop( indexX, indexY);
    prop[1]=_scgisTest->prop( indexX, indexY);
    prop[2]=_gisgpTest->prop( indexX, indexY);
    prop[3]=_scgisgpTest->prop( indexX, indexY);
    return prop;
}
 vector<double> Test::getconvAll(int ind){
	assert(_case==4);
	vector<double> conv(4);
	conv[0]=_gisTest->getconv(ind);
    conv[1]=_scgisTest->getconv(ind);
    conv[2]=_gisgpTest->getconv(ind);
    conv[3]=_scgisgpTest->getconv(ind);
    return conv;
}
vector<int> Test::getsizeconvAll(){
    assert(_case==4);
    vector<int> convsize(4);
    convsize[0]= _gisTest->getsizeconv();
    convsize[1]= _scgisTest->getsizeconv();
    convsize[2]= _gisgpTest->getsizeconv();
    convsize[3]= _scgisgpTest->getsizeconv();
    return convsize;
 }
*/
