#include "TestMI.h"

TestMI::TestMI(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, IsParameter param, int version){
	    assert(eX.columns() ==2 && eY.columns()==1);
	    vector<vector<int> > systAX(1,vector<int>(0));
	    systAX[0].push_back(0);
	    systAX[0].push_back(1);
	    vector<vector<int> > systAY(1,vector<int>(0));
	    systAY[0].push_back(0);
	    vector<vector<int> > systBX(1,vector<int>(0));
	    systBX[0].push_back(1);
	    vector<vector<int> > systBY(1,vector<int>(0));
	    systBY[0].push_back(0);
	    _cmi  = false;
        _valXUL = NULL;
        _valYUL = NULL;
	    _valX = &eX;
	    _valY = &eY;
	    _X    = &aX;
	    _Y    = &aY;
	    switch(version){
	      case 0:
	  		  _p1   = new GIS(*_valX,*_valY,*_X,*_Y,systAX,systAY, param);
	          _p2   = new GIS(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
	        break;
	      case 1:
	  		  _p1   = new SCGIS(*_valX, *_valY, *_X, *_Y, systAX, systAY, param);
	          _p2   = new SCGIS(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
	        break;
	      case 2:
	  		  _p1   = new GISgp(*_valX, *_valY, *_X, *_Y, systAX, systAY, param);
	          _p2   = new GISgp(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
	        break;
	      case 3:
	  		  _p1   = new SCGISgp(*_valX, *_valY, *_X, *_Y, systAX, systAY, param);
	          _p2   = new SCGISgp(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
	        break;
	      default:
	        cout << "default " << endl;
	  		  _p1   = new GIS(*_valX, *_valY, *_X, *_Y, systAX, systAY, param);
	          _p2   = new GIS(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
	    }

}
TestMI::TestMI(ULContainer &eX, ULContainer &eY, int version){
    assert(eX.columns() ==2 && eY.columns()==1);
    vector<vector<int> > systAX(1,vector<int>(0));
    systAX[0].push_back(0);
    systAX[0].push_back(1);
    vector<vector<int> > systAY(1,vector<int>(0));
    systAY[0].push_back(0);
    vector<vector<int> > systBX(1,vector<int>(0));
    systBX[0].push_back(1);
    vector<vector<int> > systBY(1,vector<int>(0));
    systBY[0].push_back(0);
    _cmi=true;

    _valXUL = &eX;
    _valYUL = &eY;
    _valX=NULL;
    _valY=NULL;
    _Y = NULL;
    _X=new DContainer(10,1);
    for(int i=0;i<10;i++){
    	(*_X) << i;
    }
    IsParameter param;
    param.lambdavalue    = 1.0;
    param.lambdadeltaval = 1.0;
    param.sigma          = 0.01;   //TODO find best sigma value
    param.maxit          = 10;
    param.konv           = 0.000001;
    param.time           = true;
    param.test           = true;
    param.seconds        = 10;

    switch(version){
      case 0:
  		  _p1   = new GIS(*_valXUL,*_valYUL, *_X, *_X,systAX,systAY, param);
          _p2   = new GIS(*_valXUL, *_valYUL, *_X, *_X, systBX, systBY, param);
        break;
      case 1:
  		  _p1   = new SCGIS(*_valXUL, *_valYUL, *_X, *_X, systAX, systAY, param);
          _p2   = new SCGIS(*_valXUL, *_valYUL, *_X, *_X, systBX, systBY, param);
        break;
      case 2:
  		  _p1   = new GISgp(*_valXUL, *_valYUL, *_X, *_X, systAX, systAY, param);
          _p2   = new GISgp(*_valXUL, *_valYUL, *_X, *_X, systBX, systBY, param);
        break;
      case 3:
  		  _p1   = new SCGISgp(*_valXUL, *_valYUL, *_X, *_X, systAX, systAY, param);
          _p2   = new SCGISgp(*_valXUL, *_valYUL, *_X, *_X, systBX, systBY, param);
        break;
      default:
        cout << "default " << endl;
  		  _p1   = new GIS(*_valXUL, *_valYUL, *_X, *_X, systAX, systAY, param);
          _p2   = new GIS(*_valXUL, *_valYUL, *_X, *_X, systBX, systBY, param);
    }
}

double TestMI::getMI(){
	double val;
	double p;
	if(_cmi==false){
		for(int i=0;i< pow(_X->rows(),_valX->columns() );i++ ){
	      for( int j=0; j< pow(_Y->rows(),_valY->columns() );j++){
	    	  p=_p1->propAlphX(i,j);
	    	  val+=p*_p1->propm(i)*log2(p/ (_p2->propAlphX(i,j)));
	      }
		}
		return val;
	}
	else{
		for(int i=0;i< pow(_X->rows(),_valXUL->columns() );i++ ){
	      for( int j=0; j< pow(_X->rows(),_valYUL->columns() );j++){
	    	  p=_p1->propAlphX(i,j);
	    	  val+=p*_p1->propm(i)*log2(p/ (_p2->propAlphX(i,j)));
	      }
		}
		return val;
	}

}
