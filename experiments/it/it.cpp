#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_MI.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>
#include <entropy++/sparse/state/MC_MI.h>
#include <entropy++/Matrix.h>
#include <entropy++/SparseMatrix.h>

#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"
#include "FeatureMatrix.h"
#include "GIS.h"




int main(int argc, char **argv)
{


  /*IContainer *X = new IContainer(20,1);
    Feature **FA;
  FA = new Feature*[8];
  for(int m=0; m< 8; m++){
	  FA[m]= new Feature[7];
  }
  for(int m=0;m<8; m++){
	  for(int k=0;k<8;k++){
		  FA[m][k]=new Feature(*X,*Y,(*eX)(m,0),(*eY)(k,0));
	  }
  }


  *X << 0 << 1 << 1 << 0 << 1 <<1 <<0 << 1;

  double **a;
  a = new double*[5];
  for(int m=0; m<5; m++)
	  a[m]= new double[3];

  a[0][1]=1;
  int *c= new int[3];
 int z=(*X).rows();
  cout << *X  << endl;
  cout << z << endl;
  cout << "hier" << endl;
  IContainer *T = new IContainer(8,7);
  *T << 3 << 3 << 4 << 4 << 3 << 4 << 3;
  *T << 1 << 1 << 2 << 2 << 1 << 2 << 1;
  *T << 1 << 1 << 2 << 2 << 1 << 2 << 1;
  *T << 3 << 3 << 4 << 4 << 3 << 4 << 3;
  *T << 1 << 1 << 2 << 2 << 1 << 2 << 1;
  *T << 1 << 1 << 2 << 2 << 1 << 2 << 1;
  *T << 3 << 3 << 4 << 4 << 3 << 4 << 3;
  *T << 1 << 1 << 2 << 2 << 1 << 2 << 1;
  cout << *T << endl;
    IContainer *Y = new IContainer(7,1);
    *Y << 1 << 1 << 0 << 0 << 1 << 0 << 1;

   //cout << (*X)(1,1) << endl;

	DContainer *X = new DContainer(10,1);
	*X << 1 << 2 << 3;
	DContainer *Y = new DContainer(7,1);
	*Y << 1 << 2 << 0 ;
	cout << *Y  << endl;
	DContainer *eX = new DContainer(8,1);
		*eX << 1 << 2 << 3 << 4 << 5;
	DContainer *eY = new DContainer(7,1);
		*eY << 8 << 2 << 2;
		cout << *eX  << endl;
		cout << *eY  << endl;
  Feature *F=  new Feature(1,-1);
  Feature *E= new Feature(*X,*Y,0,0);
  (*F).setlambda(1,1,2);
  double m = (*E).value();
  double k = (*F).getlambda(0,1);
  double l = (*F).value();
  cout << l << endl;
  cout << k << endl;
  cout << m << endl;
  double **a;
  a = new double*[5];
  for(int m=0; m<5; m++)
	  a[m]= new double[3];

  }*/
	vector<vector<int> > W(2,vector<int>(0));
	vector<vector<int> > ma[1][2];
	for(int i=0; i<1;i++){
		for(int j=0;j<2;j++){
			ma[i][j]=W;
		}
	}
	ma[0][1][0].push_back(2);
	//cout << ma[0][1][0][0] << endl;

  DContainer *X = new DContainer(2,1);
  		*X << 0 << 1;

  DContainer *Y = new DContainer(2,1);
  	  	*Y << 0 << 1;
  DContainer *eX = new DContainer(4,4);
  *eX << 1 << 1 << 1 << 0 << 0 << 0 << 1 << 0 << 1 << 1 << 0 << 0 <<1 << 1 << 1 << 0;
  	cout << (*eX) << endl;
  DContainer *eY = new DContainer(8,3);
  *eY << 0 << 1 << 1 << 0 << 0 << 1 << 0  << 0  << 1 << 0 << 0 << 1;
  	  	int size= (*eX).rows();
  	  	cout << size << endl;
  	 cout << (*eY) << endl;

  //FeatureMatrix *FM= new FeatureMatrix(*eX,*eY,*X,*Y,1);


  //int j= FM->getMatrixIndexX(0,0)[1];
  //int p= (*FM).getMatrixIndexdX(0,0)[0];



  /*
  double m= FM->getFeatureArrayvalue(1,1,1,0);
  cout << m << endl;
  double** observed;
  observed = new double*[2];
  for(int i=0; i<2; i++){
	  observed[i]=new double[2];
  }

	//vector observed
	for(int i=0;i<2;i++ ){
		for(int j=0; j< 2;j++){
			for(int k=0; k< FM->getMatrixIndexX(i,j).size();k++){
				observed[FM->getMatrixIndexX(i,j)[k]][FM->getMatrixIndexY(i,j)[k]]++;
			}
		}
	}

	int currmax=0;
	int curr=0;
	for(int i=0; i< 4;i++){
		for(int j=0; j< 2;j++){

			for(int featxi=0; featxi < 2; featxi++){
				for(int featyj=0; featyj < 2; featyj++){
					for(int k=0; k< FM->getMatrixIndexX(i,j).size();k++){
						if(FM->getMatrixIndexX(i,j)[k]==featxi && FM->getMatrixIndexY(i,j)[k]==featyj){
							curr++;


						}
					}
				}
			}
		if(curr> currmax) currmax=curr;
		curr=0;
		}
	}
	const int c= currmax;


	cout << c << endl;
	cout << "hier2" << endl;
	 */
	double**** observed;
	observed = new double***[2];
	for(int i=0; i<2; i++){
		observed[i]=new double**[2];
		for( int j=0;j< 2;j++){
			observed[i][j]=new double*[2];
			for(int k=0; k< 2; k++){
				observed[i][j][k]= new double[2];
				for(int l=0; l< 2;l++){
					observed[i][j][k][l]=0;
				}
			}
		}
	}
//cout << "hier " << endl;

//GIS *G = new GIS(*eX,*eY,*X,*Y,2,250,0.1);
double pq;
//pq =(*G).gis(1,2,1,0);
for(int i=0; i< 4; i++){
	for(int j=0;j<3;j++){
		//cout << (*G).gis(i,j,0,1)<< endl;
		//cout << (*G).gis(i,j,1,0)<< endl;
		//cout << (*G).gis(i,j,1,1)<< endl;
		//cout << (*G).gis(i,j,0,0)<< endl;
		//cout << endl;
	}
}
int x=rand();
srand (time(NULL));
int y= rand() % 2;
cout << y << endl;
DContainer *zX = new DContainer(2,1);
*zX << 0 << 1;
 DContainer *zY = new DContainer(2,1);
 *zY << 0 << 1;
 DContainer *zeX = new DContainer(4,5);
 for(int i=0;i< 4;i++ ){
	 for(int j=0;j<5;j++){
		 *zeX << rand() % 2;
	 }
 }
 cout << (*zeX) << endl;
 DContainer *zeY = new DContainer(4,4);
 for(int i=0;i< 4;i++ ){
	 for(int j=0;j<4;j++){
		 *zeY << rand() % 2;
	 }
 }
 FeatureMatrix *FM= new FeatureMatrix(*eX,*eY,*X,*Y,1);
 cout << (*zeY) << endl;
 GIS *Test = new GIS(*zeX,*zeY,*zX,*zY,0.1,20,0.1);

 for(int i=0; i< 5; i++){
 	for(int j=0;j<4;j++){
 		cout << (*Test).gis(i,j,0,1)<< endl;
 		cout << (*Test).gis(i,j,1,0)<< endl;
 		cout << (*Test).gis(i,j,1,1)<< endl;
 		cout << (*Test).gis(i,j,0,0)<< endl;
 		cout << endl;
 	}
 }
 //cout << (double)rand()/RAND_MAX << endl;
}
