#include "Comp.h"

Comp::Comp(GIS &exact, DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,int maxit, double konv){
	_exact= &exact;
    _valX= &eX;
    _valY= &eY;
    _X= &aX;
    _Y= &aY;
    __comptime(maxit,konv);
}
void Comp::__comptime(int maxit, double konv){
	time_t befor1;
	time_t befor2;
	time_t after1;
	time_t after2;
	befor1=time(NULL);
	_gisTest=new GIS(*_valX,*_valY,*_X,*_Y,1,maxit,konv,false);
	after1=time(NULL);
	_timediff.push_back(difftime(after1,befor1));
	befor2=time(NULL);
	_scgisTest=new SCGIS(*_valX,*_valY,*_X,*_Y,1,maxit,konv,false);
	after2=time(NULL);
	_timediff.push_back(difftime(after2,befor2));
	cout << _timediff[0] << " " << after1 << " " << befor1 << " " << endl;
	cout << _timediff[1] << endl;
	cout << "hier" << endl;
	cout <<_gisTest->getFeatureArraylambda(0,0,0,0) << endl;
	cout <<_gisTest->getFeatureArraylambda(0,0,1,0) << endl;
	cout <<_gisTest->getFeatureArraylambda(0,0,0,1) << endl;
	cout <<_gisTest->getFeatureArraylambda(0,0,1,1) << endl;
	cout << endl;
	 cout <<_scgisTest->getFeatureArraylambda(0,0,0,0) << endl;
	 cout <<_scgisTest->getFeatureArraylambda(0,0,1,0) << endl;
	 cout <<_scgisTest->getFeatureArraylambda(0,0,0,1) << endl;
	 cout <<_scgisTest->getFeatureArraylambda(0,0,1,1) << endl;
}
vector<vector<double> > Comp::getY(){
	vector<double> fill(0);
	vector<vector<double> > Y(_valY->columns(),vector<double>(pow(_Y->rows(),_valY->columns())));
	cout << Y.size() << endl;
	cout << Y[1].size() << endl;
	cout << Y[0].size() << endl;
	__fill(fill,0,Y);
	return Y;
}

void  Comp::  __fill(vector<double> fill, int i, vector<vector<double> > &Y){
	cout << "hier 1" << endl;
	cout << i << endl;
	cout << _valY->columns() << endl;
	if(_valY->columns()==i){
		static int j=0;								// nicht schoen, aber es funktioniert
		cout << "ende!!!!!!!" << endl;
		for(int l=0;l<_valY->columns();l++){
		Y[l][j]=fill[l];
		cout << Y[l][j] << " " << l << " " << j << endl;
		}
		j++;
		fill.pop_back();
		i--;
	}
	else{
		cout << " else " << endl;
		for(int x=0;x<_Y->rows();x++){
			fill.push_back(_Y->get(x,0));
			__fill(fill,i+1,Y);
			fill.pop_back();
		}

	}

}
