#include "Comp.h"

Comp::Comp(GIS &exact, DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,int maxit, double konv){
	_exact= &exact;
    _valX= &eX;
    _valY= &eY;
    _X= &aX;
    _Y= &aY;
    _sizeColValY=_valY->columns();
    _sizeColValX=_valX->columns();
    __comptime(maxit,konv);
    _alphY=__getY();
    _alphX=__getX();

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
}
void Comp:: comparison(int RowY){
		cout << "Comparison: " << endl;
		cout << endl;
		cout << "time:  GIS: " << _timediff[0] << "s SCGIS " << _timediff[1] << "s" << endl;
		cout << endl;
		vector<double> kl = KL(_alphY,RowY);
		cout << "KL-distance: GIS: " << kl[0] << " SCGIS " << kl[1] << endl;

		/*cout << "alphY :" << endl;
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
		cout << _alphX.size() << endl; */
}
vector<double> Comp:: KL(vector<vector<double> > y, int RowY ){
	vector<double> dist(2);
	double p1=0;
	double p2=0;
	double q=0;
	for(int RowX=0;RowX<  _alphX.size() ;RowX++){
		p1=_gisTest->prop(_alphX,RowX,y,RowY);
		p2=_scgisTest->prop(_alphX,RowX,y,RowY);
		q= _exact->prop(_alphX,RowX,y,RowY);
		if(fabs(p1)<0.00000001){ p1=0.000001;}
		if(fabs(p2)<0.00000001){ p2=0.000001;}
		if(fabs( q)<0.00000001){ q =0.000001;}
		dist[0]+= p1*log(p1/q);
		dist[1]+= p2*log(p2/q);
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
