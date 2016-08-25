#include "GIS.h"

GIS::GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue) {
			DContainer *valX= &eX;
			DContainer *valY= &eY;
			DContainer *X= &aX;
			DContainer *Y= &aY;
			_sizeX = (*X).rows();
			_sizeY = (*Y).rows();
			_sizeColValY= (*valY).columns();
			_sizeColValX= (*valX).columns();
			_sizeRowValX= (*valX).rows();
			_sizeRowValY= (*valY).rows();
			FeatureMatrix *FM=new FeatureMatrix(*valX,*valY,*X,*Y,lambdavalue);
			gis(*FM);
}

GIS::~GIS() {
}
double**** GIS:: __getobs(FeatureMatrix &FM){
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
				for(int k=0; k< FM.getMatrixIndexX(i,j).size();k++){
					observed[FM.getMatrixIndexX(i,j)[k]][FM.getMatrixIndexY(i,j)[k]][FM.getMatrixIndexdX(i,j)[k]][FM.getMatrixIndexdY(i,j)[k]]++;
				}
			}
		}
		return observed;

}
double** GIS::__getFeatconst(FeatureMatrix &FM){
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
							for(int k=0; k< FM.getMatrixIndexX(i,j).size();k++){
								if(FM.getMatrixIndexdX(i,j)[k]==deltxi && FM.getMatrixIndexdY(i,j)[k]==deltyj){
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
void GIS:: __getexp(FeatureMatrix &FM, double**** &expect, double*** &exponent,double** &normaliser){
	for(int Feati=0; Feati< _sizeColValX; Feati++){
		for(int Featj=0; Featj< _sizeColValY; Featj++){

			for(int xi=0; xi< _sizeRowValX; xi++){
				normaliser[Feati][Featj]=0;
				for(int yj=0; yj< _sizeRowValY; yj++){
					exponent[Feati][Featj][yj]=0;
					for(int k=0; k< FM.getMatrixIndexX(xi,yj).size();k++){
						if(FM.getMatrixIndexX(xi,yj)[k]==Feati && FM.getMatrixIndexY(xi,yj)[k]==Featj){
							exponent[Feati][Featj][yj]+= FM.getFeatureArrayvalue(Feati,Featj,xi,yj);
							if(FM.getMatrixIndexX(xi,yj)[k+1]==Feati && FM.getMatrixIndexY(xi,yj)[k+1]==Featj) k++;
						}
					}
				normaliser[Feati][Featj]=normaliser[Feati][Featj]+exp(exponent[Feati][Featj][yj]);
				}
				for(int yj=0; yj< _sizeRowValY; yj++){
					for(int k=0; k< FM.getMatrixIndexX(xi,yj).size();k++){
						if(FM.getMatrixIndexX(xi,yj)[k]==Feati && FM.getMatrixIndexY(xi,yj)[k]==Featj){
						expect[Feati][Featj][FM.getMatrixIndexdX(xi,yj)[k]][FM.getMatrixIndexdY(xi,yj)[k]]+= exp(exponent[Feati][Featj][yj])/normaliser[Feati][Featj];
						}
					}

				}
			}
		}
	}

}
void GIS:: gis(FeatureMatrix &FM){
	//observed
	double**** observ = __getobs(FM);

	//constant c for delta
	double** featconst = __getFeatconst(FM);

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
	for(int i=0; i<10;i++){
		cout << "neu" << endl;
		__getexp(FM,expected,exponent,normaliser);
		for(int Feati=0; Feati<_sizeColValX;Feati++){
			for(int Featj=0; Featj< _sizeColValY;Featj++){
				for(int lambdai=0; lambdai< _sizeX; lambdai++){
					for(int lambdaj=0; lambdaj< _sizeY; lambdaj++){
						double oldl= FM.getFeatureArraylambda(Feati,Featj,lambdai,lambdaj);
						double newl;
							if(expected[Feati][Featj][lambdai][lambdaj]!=0){
								newl= oldl + (1/featconst[Feati][Featj])*log(observ[Feati][Featj][lambdai][lambdaj]/expected[Feati][Featj][lambdai][lambdaj]);
							}
							else{
								newl=0; //
							}
							FM.setFeatureArraylambda(Feati,Featj,lambdai,lambdaj,newl);
							cout<< FM.getFeatureArraylambda(Feati,Featj,lambdai,lambdaj) << endl;
					}
				}
			}
		}
	}


}

