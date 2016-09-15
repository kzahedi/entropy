#include "Comp.h"

//vergleichswerte, gemessene X,Y und Eingabealphabete
Comp::Comp(int ColX,int RowX,int ColValY,  vector<double> lambda,DContainer &aX, DContainer &aY,int maxit, double konv,bool time,int seconds){
	_timetest=time;
    _X= &aX;
    _Y= &aY;
    _sizeColValY=ColValY;
    __getValX(ColX,RowX);
    _sizeColValX=_valX->columns();
    _exact= new GIS(ColValY,*_valX,*_X,*_Y);
    __setlambdarand(lambda);
    _alphY=__getY();
    _alphX=__getX();
    __getValY(ColValY,RowX);
    if(_timetest){
    	__comptime(maxit,konv,seconds);
    }
    else{
    	__comptime(maxit,konv);
    }
}
Comp::	~Comp(){
	_timediff.clear();
	for(int i=0;i<_alphX.size();i++){
		_alphX[i].clear();
	}
	_alphX.clear();

	for(int i=0;i<_alphY.size();i++){
		_alphY[i].clear();
	}
	_alphY.clear();
	_exact->~GIS();
	_gisTest->~GIS();
	_scgisTest->~SCGIS();
	_scgisgpTest->~SCGISgp();
	_gisgpTest->~GISgp();
}
//GIS und SCGIS ausfuehren mit Zeitmessung
void Comp::__comptime(int maxit, double konv,int seconds){
	time_t befor1;
	time_t befor2;
	time_t after1;
	time_t after2;
	time_t befor3;
	time_t after3;
	time_t befor4;
	time_t after4;
	befor1=time(NULL);
	_gisTest=new GIS(*_valX,*_valY,*_X,*_Y,1,maxit,konv,false,true,seconds);
	after1=time(NULL);
	_timediff.push_back(difftime(after1,befor1));
	befor2=time(NULL);
	_scgisTest=new SCGIS(*_valX,*_valY,*_X,*_Y,1,maxit,konv,false,true,seconds);
	after2=time(NULL);
	_timediff.push_back(difftime(after2,befor2));
	befor3=time(NULL);
	_gisgpTest=new GISgp(*_valX,*_valY,*_X,*_Y,1,1,0.01,maxit,konv,false,true,seconds);
	after3=time(NULL);
	_timediff.push_back(difftime(after3,befor3));
	befor4=time(NULL);
	_scgisgpTest= new SCGISgp(*_valX,*_valY,*_X,*_Y,1,1,0.01,maxit,konv,false,true,seconds);
	after4=time(NULL);
	_timediff.push_back(difftime(after4,befor4));
}
void Comp::__comptime(int maxit, double konv){
	time_t befor1;
	time_t befor2;
	time_t after1;
	time_t after2;
	time_t befor3;
	time_t after3;
	time_t befor4;
	time_t after4;
	befor1=time(NULL);
	_gisTest=new GIS(*_valX,*_valY,*_X,*_Y,1,maxit,konv,false,false,0);
	after1=time(NULL);
	_timediff.push_back(difftime(after1,befor1));
	befor2=time(NULL);
	_scgisTest=new SCGIS(*_valX,*_valY,*_X,*_Y,1,maxit,konv,false,false,0);
	after2=time(NULL);
	_timediff.push_back(difftime(after2,befor2));
	befor3=time(NULL);
	_gisgpTest=new GISgp(*_valX,*_valY,*_X,*_Y,1,1,0.01,maxit,konv,false,false,0);
	after3=time(NULL);
	_timediff.push_back(difftime(after3,befor3));
	befor4=time(NULL);
	_scgisgpTest= new SCGISgp(*_valX,*_valY,*_X,*_Y,1,1,0.01,maxit,konv,false,false,0);
	after4=time(NULL);
	_timediff.push_back(difftime(after4,befor4));
}
//Ausgabe
void Comp:: comparison(){
		cout << "Comparison: " << endl;
		cout << endl;
		cout << "time:  GIS: " << _timediff[0] << "s SCGIS: " << _timediff[1] << "s GIS smoothed: " << _timediff[2]<<  "s SCGIS smoothed: "  << _timediff[3] << "s" << endl;
		cout << endl;
		vector<double> kl = KL(_alphY);
		cout << "KL-distance: GIS: " << kl[0] <<  " SCGIS: " << kl[1] <<" GIS smoothed: "<< kl[2] << " SCGIS smoothed: " << kl[3] <<endl;
		if(_timetest){
		cout<< "Iterations: GIS: " << _gisTest->getIterations() << " SCGIS: " << _scgisTest->getIterations()<< " GIS smoothed: " << _gisgpTest->getIterations() << " SCGIS smoothed: " << _scgisgpTest->getIterations() <<   endl;
		}
		cout << "lambda: " << endl;
		cout << "vergleichswerte" << endl;
		 cout << _exact->getFeatureArraylambda(0,0,0,0) <<endl;
		 cout << _exact->getFeatureArraylambda(0,0,1,0) <<endl;
		 cout << _exact->getFeatureArraylambda(0,0,0,1) <<endl;
		 cout << _exact->getFeatureArraylambda(0,0,1,1) <<endl;
		 cout << endl;
		 cout << _exact->getFeatureArraylambda(1,0,0,0) <<endl;
		 cout << _exact->getFeatureArraylambda(1,0,1,0) <<endl;
		 cout << _exact->getFeatureArraylambda(1,0,0,1) <<endl;
		 cout << _exact->getFeatureArraylambda(1,0,1,1) <<endl;
		 cout << endl;
		 cout << _exact->getFeatureArraylambda(0,1,0,0) <<endl;
		 cout << _exact->getFeatureArraylambda(0,1,1,0) <<endl;
		 cout << _exact->getFeatureArraylambda(0,1,0,1) <<endl;
		 cout << _exact->getFeatureArraylambda(0,1,1,1) <<endl;
		 cout << endl;
		 cout << _exact->getFeatureArraylambda(1,1,0,0) <<endl;
		 cout << _exact->getFeatureArraylambda(1,1,1,0) <<endl;
		 cout << _exact->getFeatureArraylambda(1,1,0,1) <<endl;
		 cout << _exact->getFeatureArraylambda(1,1,1,1) <<endl;

		cout << "GIS" << endl;
		 cout << _gisTest->getFeatureArraylambda(0,0,0,0) <<endl;
		 cout << _gisTest->getFeatureArraylambda(0,0,1,0) <<endl;
		 cout << _gisTest->getFeatureArraylambda(0,0,0,1) <<endl;
		 cout << _gisTest->getFeatureArraylambda(0,0,1,1) <<endl;
		 cout << endl;
		 cout << _gisTest->getFeatureArraylambda(1,0,0,0) <<endl;
		 cout << _gisTest->getFeatureArraylambda(1,0,1,0) <<endl;
		 cout << _gisTest->getFeatureArraylambda(1,0,0,1) <<endl;
		 cout << _gisTest->getFeatureArraylambda(1,0,1,1) <<endl;
		 cout << endl;
		 cout << _gisTest->getFeatureArraylambda(0,1,0,0) <<endl;
		 cout << _gisTest->getFeatureArraylambda(0,1,1,0) <<endl;
		 cout << _gisTest->getFeatureArraylambda(0,1,0,1) <<endl;
		 cout << _gisTest->getFeatureArraylambda(0,1,1,1) <<endl;
		 cout << endl;
		 cout << _gisTest->getFeatureArraylambda(1,1,0,0) <<endl;
		 cout << _gisTest->getFeatureArraylambda(1,1,1,0) <<endl;
		 cout << _gisTest->getFeatureArraylambda(1,1,0,1) <<endl;
		 cout << _gisTest->getFeatureArraylambda(1,1,1,1) <<endl;
		cout << "GIS smoothed " << endl;
		 cout << _gisgpTest->getFeatureArraylambda(0,0,0,0) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(0,0,1,0) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(0,0,0,1) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(0,0,1,1) <<endl;
		 cout << endl;
		 cout << _gisgpTest->getFeatureArraylambda(1,0,0,0) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(1,0,1,0) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(1,0,0,1) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(1,0,1,1) <<endl;
		 cout << endl;
		 cout << _gisgpTest->getFeatureArraylambda(0,1,0,0) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(0,1,1,0) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(0,1,0,1) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(0,1,1,1) <<endl;
		 cout << endl;
		 cout << _gisgpTest->getFeatureArraylambda(1,1,0,0) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(1,1,1,0) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(1,1,0,1) <<endl;
		 cout << _gisgpTest->getFeatureArraylambda(1,1,1,1) <<endl;

		cout << "alphY :" << endl;
		for(int i=0;i<_alphY.size();i++){
			for(int j=0;j<_alphY[i].size();j++){
				cout << _alphY[i][j] << " " ;
			}
			cout << endl;
		}
		cout << "alphX :" << endl;
		for(int i=0;i<_alphX.size();i++){
			for(int j=0;j<_alphX[i].size();j++){
				cout << _alphX[i][j] << " " ;
			}
			cout << endl;
		}
		cout << _alphX.size() << endl;

}
//Abstand
vector<double> Comp:: KL(vector<vector<double> > y ){
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
	for(int RowX=0;RowX<  _alphX.size() ;RowX++){
		pm1=_gisTest->propm(_alphX,RowX,y);
		pm2=_scgisTest->propm(_alphX,RowX,y);
		pm3=_gisgpTest->propm(_alphX,RowX,y);
		pm4=_scgisgpTest->propm(_alphX,RowX,y);
		for(int sizeY=0;sizeY<y.size();sizeY++){
			p1=_gisTest->propm(_alphX,RowX,y,sizeY);
			p2=_scgisTest->propm(_alphX,RowX,y,sizeY);
			p3=_gisgpTest->propm(_alphX,RowX,y,sizeY);
			p4=_scgisgpTest->propm(_alphX,RowX,y,sizeY);
			q= _exact->propm(_alphX,RowX,y,sizeY);
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

vector<vector<double> > Comp::__getY(){
	vector<double> fill(0);
	vector<vector<double> > Y(pow(_Y->rows(),_sizeColValY),vector<double>(_sizeColValY));
	__fill(fill,0,Y);
	return Y;
}
vector<vector<double> > Comp::__getX(){
	vector<double> fillx(0);
	vector<vector<double> > X(pow(_X->rows(),_sizeColValX),vector<double>(_sizeColValX));
	__fillx(fillx,0,X);
	return X;
}
// nicht schoen, aber es funktioniert
void  Comp::  __fill(vector<double> fill, int i, vector<vector<double> > &Y){
	if(_sizeColValY==i){
		static int j=0;
		for(int l=0;l<_sizeColValY;l++){
			Y[j][l]=fill[l];
		}
		j++;
		fill.pop_back();
		i--;
	}
	else{
		for(int x=0;x<_Y->rows();x++){
			fill.push_back(_Y->get(x,0));
			__fill(fill,i+1,Y);
			fill.pop_back();
		}
	}
}
void  Comp:: __fillx(vector<double> fillx, int i, vector<vector<double> > &X){
	if(_sizeColValX==i){
		static int j=0;
		for(int l=0;l<_sizeColValX;l++){
			X[j][l]=fillx[l];
		}
		j++;
		fillx.pop_back();
		i--;
	}
	else{
		for(int x=0;x<_X->rows();x++){
			fillx.push_back(_X->get(x,0));
			__fillx(fillx,i+1,X);
			fillx.pop_back();
		}
	}
}
//DContainer mit Indizes fuer das lambda und den Wert
void Comp::	__setlambda(IContainer &indizes, DContainer &values ){
	assert(indizes.columns() ==4);
	assert(indizes.rows() == values.rows());
	for(int i=0;i< indizes.rows();i++){
		_exact->setFeatureArraylambda(indizes.get(i,0),indizes.get(i,1),indizes.get(i,2),indizes.get(i,3),values.get(i,4));
	}
}
void Comp:: __setlambdarand(vector<double> lambda){
	srand(time(NULL));
	for(int Feati=0;Feati<_sizeColValX;Feati++){
		for(int Featj=0; Featj<_sizeColValY;Featj++){
			for(int delti=0;delti<_X->rows();delti++){
				for(int deltj=0; deltj<_Y->rows();deltj++){
					double r = rand() % lambda.size();
					_exact->setFeatureArraylambda(Feati,Featj,delti,deltj,lambda[r]);
				}
			}
		}
	}
}
void Comp::__getValX(int ColX,int RowX){
	_valX = new DContainer(RowX,ColX);
	for(int i=0;i< RowX;i++ ){
		for(int j=0;j<ColX;j++){
			double z = rand() % _X->rows();
			*_valX << (*_X)(z,0) ;
		}
	}
}
void	Comp::__getValY(int ColY,int RowX){
	srand(time(NULL));
	 double** prop;
	 prop=new double*[RowX];
	 for(int i=0;i<RowX;i++){
		 prop[i]=new double[_alphY.size()];
		 for(int j=0;j<_alphY.size();j++){
			 prop[i][j]=0;
		 }
	 }
	for(int i=0;i<RowX;i++){
		for(int propi=0;propi<_alphY.size();propi++){
			double l=_exact->prop(i,_alphY,propi);
				prop[i][propi]=_exact->prop(i,_alphY,propi);
			}
	}
	_valY = new DContainer(RowX,ColY);
	for(int i=0;i<RowX;i++){
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
}
