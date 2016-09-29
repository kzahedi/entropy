#include "TestMI.h"

TestMI::TestMI(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<vector<int> > systAX, vector<vector<int> > systAY, vector<vector<int> > systBX, vector<vector<int> > systBY, IsParameter param, int version){
	    _valX = &eX;
	    _valY = &eY;
	    _X    = &aX;
	    _Y    = &aY;
	    _exact =NULL;
	    switch(version){
	      case 0:
	  		  _p1   = new GIS(*_valX, *_valY, *_X, *_Y, systAX, systAY, param);
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
TestMI::TestMI(int colX,int colValY, int rowX,DContainer &aX, DContainer &aY,vector<double> lambda, vector<vector<int> > systAX, vector<vector<int> > systAY, vector<vector<int> > systBX, vector<vector<int> > systBY,vector<vector<int> > systCX, vector<vector<int> > systCY, IsParameter param){
	  _X           = &aX;
	  _Y           = &aY;
	  _valX = new DContainer(rowX,colX);
	  for(int i=0;i < rowX; i++)
	  {
	    for(int j=0;j < colX; j++)
	    {
	      double z = rand() % _X->rows(); // random indices for the alphabet
	      *_valX << (*_X)(z,0);
	    }
	  }
	   _exact    = new IT(colValY,*_valX,*_X,*_Y,systAX,systAY);
	  srand(time(NULL));
	  for(int feat=0;feat<systAX.size();feat++)
	  {
	    for(int delti=0;delti<pow(_X->rows(),systAX[feat].size());delti++)
	    {
	      for(int deltj=0; deltj<pow(_Y->rows(),systAY[feat].size());deltj++)
	      {
	        double r = rand() % lambda.size();
	        _exact->setFeatureArraylambda(feat,delti,deltj,lambda[r]);
	      }
	    }
	  }
	  srand(time(NULL));
	  _valY = new DContainer(rowX,colValY);
	  for(int i=0;i<rowX;i++)
	  {
	    double z= (double) rand()/RAND_MAX;
	    double s=0;
	    int ind=0;
	    for(int j=0;j< pow(_Y->rows(), colValY) && s<z;j++)
	    {
	      s+= _exact->prop(i,j);;
	      ind=j;
	    }
	    for(int k=0;k< colValY;k++)
	    {
	      (*_valY) << _exact->index(ind,false, colValY)[k] ;
	    }
	  }
	   _p1   = new GIS(*_valX, *_valY, *_X, *_Y, systCX, systCY, param);
       _p2   = new GIS(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
       cout << "GIS: "<<  getMI()  << endl;
	   _p1   = new SCGIS(*_valX, *_valY, *_X, *_Y, systCX, systCY, param);
       _p2   = new SCGIS(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
       cout << "SCGIS: " << getMI()  <<endl;
	   _p1   = new GISgp(*_valX, *_valY, *_X, *_Y, systCX, systCY, param);
       _p2   = new GISgp(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
       cout << "GISgp: " << getMI() <<endl;
	   _p1   = new SCGISgp(*_valX, *_valY, *_X, *_Y, systCX, systCY, param);
	   _p2   = new SCGISgp(*_valX, *_valY, *_X, *_Y, systBX, systBY, param);
	   cout << "SCGISgp: " << getMI()  <<endl;

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
