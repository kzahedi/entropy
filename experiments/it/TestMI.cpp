#include "TestMI.h"

TestMI::TestMI(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<double> lambda, vector<vector<int> > systAX, vector<vector<int> > systAY, vector<vector<int> > systBX, vector<vector<int> > systBY, IsParameter param, int version){
	    _valX = &eX;
	    _valY = &eY;
	    _X    = &aX;
	    _Y    = &aY;
	    cout << " Konstruktor" << endl;
		_p1   = new GIS(*_valX, *_valY, *_X, *_Y, systAX, systAY, param);
		cout << " nach p1 " << endl;
        _p2   = new GIS(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
}

double TestMI::getMI(){
	double val;
	double p;
	for(int i=0;i< pow(_X->rows(),_valX->columns() );i++ ){
      for( int j=0; j< pow(_Y->rows(),_valY->columns() );j++){
    	  p=_p1->prop(i,j);
    	  val+=p*log(p/ (_p2->prop(i,j)));
      }
	}
	return val;
}
