#include "GIS.h"

GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue) {
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
			__gis();
}

double GIS::gis(int Feati,int Featj,double ValX,double ValY){
	double norm=0;
	double exponent= exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,ValY) );
	//cout << exponent<< endl;
	for(int yi=0;yi<_sizeY;yi++){
		norm+= exp( (*_FM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
	}

	return exponent/norm;
}
GIS::~GIS() {
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
		for(int i=0;i<_sizeColValX;i++ ){
			for(int j=0; j< _sizeColValY;j++){
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
				for(int j=0; j< _sizeRowValY;j++){

					for(int deltxi=0; deltxi < _sizeX; deltxi++){
						for(int deltyj=0; deltyj < _sizeY; deltyj++){
							for(int k=0; k< (*_FM).getMatrixIndexX(i,j).size();k++){
								if((*_FM).getMatrixIndexdX(i,j)[k]==deltxi && (*_FM).getMatrixIndexdY(i,j)[k]==deltyj){
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
	for(int Feati=0; Feati< _sizeColValX; Feati++){
		for(int Featj=0; Featj< _sizeColValY; Featj++){

			for(int xi=0; xi< _sizeRowValX; xi++){
				normaliser[Feati][Featj]=0;
				for(int yj=0; yj< _sizeRowValY; yj++){
					exponent[Feati][Featj][yj]=0;
					for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++){
						if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj){
							double m =(*_FM).getFeatureArrayvalue(Feati,Featj,xi,yj);
							exponent[Feati][Featj][yj]+= (*_FM).getFeatureArrayvalue(Feati,Featj,xi,yj);
							if((*_FM).getMatrixIndexX(xi,yj)[k+1]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k+1]==Featj){
								k++;

							}
						}
					}
				normaliser[Feati][Featj]+=exp(exponent[Feati][Featj][yj]);
				}
				for(int yj=0; yj< _sizeRowValY; yj++){
					for(int k=0; k< (*_FM).getMatrixIndexX(xi,yj).size();k++){
						if((*_FM).getMatrixIndexX(xi,yj)[k]==Feati && (*_FM).getMatrixIndexY(xi,yj)[k]==Featj){
						expect[Feati][Featj][(*_FM).getMatrixIndexdX(xi,yj)[k]][(*_FM).getMatrixIndexdY(xi,yj)[k]]+= exp(exponent[Feati][Featj][yj])/normaliser[Feati][Featj];
						}
					}
				}
			}
		}
	}
}
void GIS:: __gis(){
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
				for(int k=0;k< _sizeRowValY;k++){
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
	for(int i=0; i<50;i++){
		__getexp(expected,exponent,normaliser);
		for(int Feati=0; Feati<_sizeColValX;Feati++){
			for(int Featj=0; Featj< _sizeColValY;Featj++){
				for(int lambdai=0; lambdai< _sizeX; lambdai++){
					for(int lambdaj=0; lambdaj< _sizeY; lambdaj++){
						double oldl= (*_FM).getFeatureArraylambda(Feati,Featj,lambdai,lambdaj);
						double newl;
							if(expected[Feati][Featj][lambdai][lambdaj]!=0 && observ[Feati][Featj][lambdai][lambdaj]!=0 ){
								newl= oldl + (1/featconst[Feati][Featj])*log(observ[Feati][Featj][lambdai][lambdaj]/expected[Feati][Featj][lambdai][lambdaj]);
							}
							else{
								newl=0; //
							}
							double m =(*_FM).getFeatureArrayvalue(Feati,Featj,1,1);
							(*_FM).setFeatureArraylambda(Feati,Featj,lambdai,lambdaj,newl);
					}
				}
			}
		}
	}

}

