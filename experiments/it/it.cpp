#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_MI.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>
#include <entropy++/sparse/state/MC_MI.h>
#include <entropy++/Matrix.h>
#include <entropy++/SparseMatrix.h>

#include <time.h>
#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"
#include "FeatureMatrix.h"
#include "GIS.h"
#include "InstanceMatrix.h"
#include "SCGIS.h"
#include "Comp.h"



int main(int argc, char **argv){
	srand(time(NULL));
		 int n=100000;
		 DContainer *eX = new DContainer(n,1);
		 for(int i=0;i< n;i++ ){
			 for(int j=0;j<1;j++){
				 *eX << rand() % 2;
			 }
		 }
		 DContainer *zX = new DContainer(2,1);
		 *zX << 0 << 1;
		  DContainer *zY = new DContainer(2,1);
		 *zY << 0 << 1;

		 GIS *Test = new GIS(1,*eX,*zX,*zY);
		 Test->setFeatureArraylambda(0,0,1,0,2);
		 Test->setFeatureArraylambda(0,0,1,1,0);
		 Test->setFeatureArraylambda(0,0,0,0,1);
		 Test->setFeatureArraylambda(0,0,0,1,3);

		 double** prop;
		 prop=new double*[n];
		 for(int i=0;i<n;i++){
			 prop[i]=new double[2];
			 for(int j=0;j<2;j++){
				 prop[i][j]=0;
			 }
		 }

		 for(int i=0;i<n;i++){
			for(int propi=0;propi<2;propi++){
					prop[i][propi]=Test->prop(0,0,(*eX)(i,0),propi);
				 }
			 }
		 DContainer *esY=new DContainer(n,1);

		 for(int i=0;i<n;i++){
			 double z=(double)rand()/RAND_MAX;
			 double s=0;
			 int ind=0;
			 for(int j=0;j<2 && s<z;j++){
					s+=prop[i][j];
					ind=j;
				 }
			 (*esY) << ind;
			 }
		 //GIS &exact, DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,int maxit, double konv
		 Comp *test = new Comp(*Test,*eX,*esY,*zX,*zY,2000,0.000001);
/*srand(time(NULL));
int n=1000;
DContainer *eX = new DContainer(n,1);
for(int i=0;i< n;i++ ){
	 for(int j=0;j<1;j++){
		 *eX << rand() % 2;
	 }
}
DContainer *zX = new DContainer(2,1);
*zX << 0 << 1;
 DContainer *zY = new DContainer(2,1);
*zY << 0 << 1;

GIS *Test = new GIS(1,*eX,*zX,*zY);
Test->setFeatureArraylambda(0,0,1,0,4);
Test->setFeatureArraylambda(0,0,1,1,0);
Test->setFeatureArraylambda(0,0,0,0,3);
Test->setFeatureArraylambda(0,0,0,1,2);

double** prop;
prop=new double*[n];
for(int i=0;i<n;i++){
	 prop[i]=new double[2];
	 for(int j=0;j<2;j++){
		 prop[i][j]=0;
	 }
}

for(int i=0;i<n;i++){
	for(int propi=0;propi<2;propi++){
			prop[i][propi]=Test->gis(0,0,(*eX)(i,0),propi);
		 }
	 }
DContainer *esY=new DContainer(n,1);

for(int i=0;i<n;i++){
	 double z=(double)rand()/RAND_MAX;
	 double s=0;
	 int ind=0;
	 for(int j=0;j<2 && s<z;j++){
			s+=prop[i][j];
			ind=j;
		 }
	 (*esY) << ind;
	 }
GIS *zTest = new GIS(*eX,*esY,*zX,*zY,1,5000,0.01,true);
SCGIS *bTest = new SCGIS(*eX,*esY,*zX,*zY,5000,0.01,1,true);
//cout << "exakt:"  << endl;
//cout <<Test->gis(0,0,0,0) << " " << Test->getFeatureArraylambda(0,0,0,0) << endl;
//cout <<Test->gis(0,0,1,0) <<  " " << Test->getFeatureArraylambda(0,0,1,0) <<endl;
//cout <<Test->gis(0,0,0,1) << " " << Test->getFeatureArraylambda(0,0,0,1) << endl;
//cout <<Test->gis(0,0,1,1) <<  " " << Test->getFeatureArraylambda(0,0,1,1) <<endl;
cout << "gis: " << endl;
cout <<zTest->gis(0,0,0,0)-Test->gis(0,0,0,0) <<  " " << zTest->getFeatureArraylambda(0,0,0,0) <<endl;
cout <<zTest->gis(0,0,1,0)-Test->gis(0,0,1,0) <<  " " << zTest->getFeatureArraylambda(0,0,1,0) <<endl;
cout <<zTest->gis(0,0,0,1)-Test->gis(0,0,0,1)  <<  " " << zTest->getFeatureArraylambda(0,0,0,1) <<endl;
cout <<zTest->gis(0,0,1,1)-Test->gis(0,0,1,1) <<  " " << zTest->getFeatureArraylambda(0,0,1,1) <<endl;
cout << "scgis " << endl;
cout <<bTest->scgis(0,0,0,0)-Test->gis(0,0,0,0) <<  " " << bTest->getFeatureArraylambda(0,0,0,0) <<endl;
cout <<bTest->scgis(0,0,1,0)-Test->gis(0,0,1,0) <<  " " << bTest->getFeatureArraylambda(0,0,1,0) <<endl;
cout <<bTest->scgis(0,0,0,1)-Test->gis(0,0,0,1)  <<  " " << bTest->getFeatureArraylambda(0,0,0,1) <<endl;
cout <<bTest->scgis(0,0,1,1)-Test->gis(0,0,1,1) <<  " " << bTest->getFeatureArraylambda(0,0,1,1) <<endl;

*/
}

/*
	vector<vector<int> > W(2,vector<int>(0));
	vector<vector<int> > ma[1][2];
	for(int i=0; i<1;i++){
		for(int j=0;j<2;j++){
			ma[i][j]=W;
		}
	}
	ma[0][1][0].push_back(2);

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

int y= rand() % 2;
//cout << y << endl;
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
 //cout << (*zeX) << endl;
 DContainer *zeY = new DContainer(4,4);
 for(int i=0;i< 4;i++ ){
	 for(int j=0;j<4;j++){
		 *zeY << rand() % 2;
	 }
 }
 srand(time(NULL));
 int n=10000;
 DContainer *eX = new DContainer(n,2);
 for(int i=0;i< n;i++ ){
	for(int j=0;j<2;j++){
		*eX << rand() % 2;
	}
 }


 GIS *Test = new GIS(2,*eX,*zX,*zY);

	 Test->setFeatureArraylambda(0,0,1,0,4);
	 Test->setFeatureArraylambda(0,0,1,1,0);
	 Test->setFeatureArraylambda(0,0,0,0,2);
	 Test->setFeatureArraylambda(0,0,0,1,5);

	 Test->setFeatureArraylambda(1,0,1,0,2);
	 Test->setFeatureArraylambda(1,0,1,1,3);
	 Test->setFeatureArraylambda(1,0,0,0,0);
	 Test->setFeatureArraylambda(1,0,0,1,6);

	 Test->setFeatureArraylambda(0,1,1,0,0);
	 Test->setFeatureArraylambda(0,1,1,1,1);
	 Test->setFeatureArraylambda(0,1,0,0,3);
	 Test->setFeatureArraylambda(0,1,0,1,5);

	 Test->setFeatureArraylambda(1,1,1,0,1);
	 Test->setFeatureArraylambda(1,1,1,1,3);
	 Test->setFeatureArraylambda(1,1,0,0,2);
	 Test->setFeatureArraylambda(1,1,0,1,1);

	 double** prop;
	 prop=new double*[n];
	 for(int i=0;i<n;i++){
		 prop[i]=new double[4];
		 for(int j=0;j<4;j++){
			 prop[i][j]=0;
		 }
	 }

	vector<vector<double> > val(4,vector<double>(2));
		 val[0][0]=0;
		 val[0][1]=0;
		 val[1][0]=1;
		 val[1][1]=0;
		 val[2][0]=0;
		 val[2][1]=1;
		 val[3][1]=1;
		 val[3][1]=1;
	for(int i=0;i<n;i++){
		for(int propi=0;propi<4;propi++){
				prop[i][propi]=Test->gis(i,val,propi);
			}
		 }
	 DContainer *esY=new DContainer(n,2);

	 for(int i=0;i<n;i++){
		 double z=(double)rand()/RAND_MAX;
		 double s=0;
		 int ind=0;
		 for(int j=0;j<4 && s<z;j++){
				s+=prop[i][j];
				ind=j;
		 }
		 if(ind==0) (*esY) << 0 << 0;
		 if(ind==1) (*esY) << 1 << 0;
		 if(ind==2) (*esY) << 0 << 1;
		 if(ind==3) (*esY) << 1 << 1;

	 }

	 GIS *zTest = new GIS(*eX,*esY,*zX,*zY,1,1000000,0.001);

		 cout << endl;
		 cout <<zTest->gis(0,0,0,0)-Test->gis(0,0,0,0)  << endl;
		 cout <<zTest->gis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
		 cout <<zTest->gis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
		 cout <<zTest->gis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
		 cout << endl;
		 cout <<zTest->gis(1,0,0,0)-Test->gis(1,0,0,0) << endl;
		 cout <<zTest->gis(1,0,1,0)-Test->gis(1,0,1,0) << endl;
		 cout <<zTest->gis(1,0,0,1)-Test->gis(1,0,0,1) << endl;
		 cout <<zTest->gis(1,0,1,1)-Test->gis(1,0,1,1) << endl;
		 cout << endl;
		 cout <<zTest->gis(0,1,0,0)-Test->gis(0,1,0,0) << endl;
		 cout <<zTest->gis(0,1,1,0)-Test->gis(0,1,1,0) << endl;
		 cout <<zTest->gis(0,1,0,1)-Test->gis(0,1,0,1) << endl;
		 cout <<zTest->gis(0,1,1,1)-Test->gis(0,1,1,1) << endl;
		 cout << endl;
		 cout <<zTest->gis(1,1,0,0)-Test->gis(1,1,0,0) << endl;
		 cout <<zTest->gis(1,1,1,0)-Test->gis(1,1,1,0) << endl;
		 cout <<zTest->gis(1,1,0,1)-Test->gis(1,1,0,1) << endl;
		 cout <<zTest->gis(1,1,1,1)-Test->gis(1,1,1,1) << endl;

*/




