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

	    _valX = &eX;
	    _valY = &eY;
	    _X    = &aX;
	    _Y    = &aY;
	    _exact =NULL;
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

double TestMI::getMI(){
	double val;
	double p;
	for(int i=0;i< pow(_X->rows(),_valX->columns() );i++ ){
      for( int j=0; j< pow(_Y->rows(),_valY->columns() );j++){
    	  p=_p1->propAlphX(i,j);
    	  val+=p*_p1->propm(i)*log2(p/ (_p2->propAlphX(i,j)));
      }
	}
	return val;
}
