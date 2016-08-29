#include "GIS.h"

GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv) {
			_valX= &eX;
			_valY= &eY;
			_X= &aX;
			_Y= &aY;
			_sizeX = (*_X).rows();
			_sizeY = (*_Y).rows();
			_sizeColValY= (*_valY).columns();
			_sizeColValX= (*_valX).columns();
			_sizeRowValX= (*_valX).rows();
			_sizeRowValY= (*_valY).rows();
			_FM=new FeatureMatrix(*_valX,*_valY,*_X,*_Y,lambdavalue);
			__gis(maxit, konv);
}
GIS::GIS(int sizeaX, int sizeaY, int sizeRowX, int sizeRowY, int sizeColX, int sizeColY ){
		_sizeX=sizeaX;
		_sizeY=sizeaY;
		_sizeColValX=sizeColX;
		_sizeColValY=sizeColY;
		_sizeRowValX=sizeRowX;
		_sizeRowValY=sizeRowY;
		_valX= new DContainer(_sizeRowValX,_sizeColValX);
		_valY= new DContainer(_sizeRowValY,_sizeRowValX);
		_X=new DContainer(_sizeX,0);
		_Y=new DContainer(_sizeY,0);
		_FM=new FeatureMatrix(*_valX,*_valY,*_X,*_Y,0);
}
double GIS::gis(int Feati,int Featj,double ValX,double ValY){
	double norm=0;
	double exponent= exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,ValY) );
	for(int yi=0;yi<_sizeY;yi++){
		norm+= exp( (*_FM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
	}

	return exponent/norm;
}
void GIS::setFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY,double valuelambda){
	//
}
double GIS::getFeatureArraylambda(int Feati, int Featj, int ilambdaX, int ilambdaY){
	//
}
double**** GIS:: __getobs(){
	double**** observed;
		observed = new double***[_sizeColValX];
		for(int i=0; i<_sizeColValX; i++){
			observed[i]=new double**[_sizeColValY];
			for( int j=0;j< _sizeColValY;j++){
				observed[i][j]=new double*[_sizeX];
				for(int k=0; k< _sizeX; k++){
					observed[i][j][k]= new double[_sizeY];
					for(int l=0; l< _sizeY;l++){
						observed[i][j][k][l]=0;
					}
				}
			}
		}

		//vector observed
		for(int i=0;i<_sizeRowValX;i++ ){
			for(int j=0; j< _sizeY;j++){
				for(int k=0; k< (*_FM).getMatrixIndexX(i,j).size();k++){
					observed[(*_FM).getMatrixIndexX(i,j)[k]][(*_FM).getMatrixIndexY(i,j)[k]][(*_FM).getMatrixIndexdX(i,j)[k]][(*_FM).getMatrixIndexdY(i,j)[k]]++;
				}
			}
		}
		return observed;


}
double** GIS::__getFeatconst(){
	double** Featconst;
	Featconst = new double*[_sizeColValX];
	for(int i=0; i< _sizeColValX;i++){
		Featconst[i]=new double[_sizeColValY];
		for(int j=0; j< _sizeColValY; j++){
			Featconst[i][j]=0;
		}
	}
	int curr=0;
	for(int delti=0; delti< _sizeColValX; delti++){
		for(int deltj=0; deltj< _sizeColValY; deltj++){

			for(int i=0; i< _sizeRowValX;i++){
				for(int j=0; j< _sizeY;j++){

					for(int deltxi=0; deltxi < _sizeX; deltxi++){
						for(int deltyj=0; deltyj < _sizeY; deltyj++){
							for(int k=0; k< (*_FM).getMatrixIndexX(i,j).size();k++){
								if((*_FM).getMatrixIndexX(i,j)[k]==delti && (*_FM).getMatrixIndexY(i,j)[k]==deltj
									&&(*_FM).getMatrixIndexdX(i,j)[k]==deltxi && (*_FM).getMatrixIndexdY(i,j)[k]==deltyj){
									curr++;
								}
							}
						}
					}
					if(curr> Featconst[delti][deltj]) Featconst[delti][deltj]=curr;
					curr=0;

				}
			}
		}
	}
	return Featconst;

}
void GIS:: __getexp(double**** &expect, double*** &exponent,double** &normaliser){

	for(int i=0; i<_sizeColValX; i++){
		for( int j=0;j< _sizeColValY;j++){
			for(int k=0; k< _sizeX; k++){
				for(int l=0; l< _sizeY;l++){
					expect[i][j][k][l]=0;
				}
			}
		}
	}
	for(int Feati=0; Feati< _sizeColValX; Feati++){
		for(int Featj=0; Featj< _sizeColValY; Featj++){

			for(int xi=0; xi< _sizeRowValX; xi++){
				normaliser[Feati][Featj]=0;
				for(int yj=0; yj< _sizeY; yj++){
					exponent[Feati][Featj][yj]=0;
					for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++){
						if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj){
							exponent[Feati][Featj][yj]+= (*_FM).getFeatureArraylambda(Feati, Featj,(*_FM).getMatrixIndexdX(xi,yj)[k], (*_FM).getMatrixIndexdY(xi,yj)[k]);
						}
					}
				//normaliser[Feati][Featj]+=exp(exponent[Feati][Featj][yj]);
				}
				//cout << normaliser[0][0] << endl;
				for(int yj=0; yj< _sizeY; yj++){
					for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++){
						if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj){
						expect[Feati][Featj][(*_FM).getMatrixIndexdX(xi,yj)[k]][(*_FM).getMatrixIndexdY(xi,yj)[k]]+=(*_FM).getFeatureArraydelta(Feati,Featj,(*_FM).getMatrixIndexdX(xi,yj)[k],(*_FM).getMatrixIndexdY(xi,yj)[k],xi,yj)*exp(exponent[Feati][Featj][yj]);
						// /normaliser[Feati][Featj]
						}
					}
				}
			}
		}
	}
}
void GIS:: __gis(int maxit, double konv){
	//observed
	double**** observ = __getobs();

	//constant c for delta
	double** featconst = __getFeatconst();

	//array for exp
	double**** expected;
		expected = new double***[_sizeColValX];
		for(int i=0; i<_sizeColValX; i++){
			expected[i]=new double**[_sizeColValY];
			for( int j=0;j< _sizeColValY;j++){
				expected[i][j]=new double*[_sizeX];
				for(int k=0; k< _sizeX; k++){
					expected[i][j][k]= new double[_sizeY];
					for(int l=0; l< _sizeY;l++){
						expected[i][j][k][l]=0;
					}
				}
			}
		}
	double*** exponent;
		exponent = new double**[_sizeColValX];
		for(int i=0; i<_sizeColValX; i++){
			exponent[i]=new double*[_sizeColValY];
			for(int j=0;j<_sizeColValY;j++){
				exponent[i][j]=new double[_sizeRowValY];
				for(int k=0;k< _sizeY;k++){
					exponent[i][j][k]=0;
				}
			}
		}
	double** normaliser;
		normaliser= new double*[_sizeColValX];
		for(int i=0; i< _sizeColValX; i++){
			normaliser[i]=new double[_sizeColValY];
			for(int j=0; j<_sizeColValY;j++){
				normaliser[i][j]=0;
			}
		}
	int i=0;
	double l=1;
	while(i<maxit && (l>=konv || l<=-konv)){
		l=0;
		__getexp(expected,exponent,normaliser);
		for(int Feati=0; Feati<_sizeColValX;Feati++){
			for(int Featj=0; Featj< _sizeColValY;Featj++){
				for(int lambdai=0; lambdai< _sizeX; lambdai++){
					for(int lambdaj=0; lambdaj< _sizeY; lambdaj++){
						double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,lambdai,lambdaj);
						double newl;
							if(expected[Feati][Featj][lambdai][lambdaj]!=0 && observ[Feati][Featj][lambdai][lambdaj]!=0 ){
								double p=(observ[Feati][Featj][lambdai][lambdaj]/expected[Feati][Featj][lambdai][lambdaj]);
								newl= oldl + (1/featconst[Feati][Featj])*log(p);
								l+=(observ[Feati][Featj][lambdai][lambdaj]-expected[Feati][Featj][lambdai][lambdaj]);
							}
							else{
								newl=0;
							}
							(*_FM).setFeatureArraylambda(Feati,Featj,lambdai,lambdaj,newl);
					}
				}
			}
		}
		i++;
	}

}

