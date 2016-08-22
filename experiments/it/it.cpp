#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_MI.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>
#include <entropy++/sparse/state/MC_MI.h>

#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"
#include "FeatureMatrix.h"


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

  DContainer *X = new DContainer(10,1);
  	*X << 1 << 1 << 1 << 1;
  DContainer *Y = new DContainer(7,1);
    *Y << 1 << 1 << 1 ;
  DContainer *eX = new DContainer(8,1);
  				*eX << 1 << 1 << 1 << 1 << 1;
  DContainer *eY = new DContainer(7,1);
  				*eY << 1 << 1 << 1;
  double a=(*X)(1,0);
  cout << a << endl;
  Feature **FA;
  FA = new Feature*[8];
for(int m=0; m< 8; m++){
  FA[m]= new Feature[7];
}
for(int m=0;m<8; m++){
  for(int k=0;k<7;k++){
	  Feature *K=new Feature(*X,*Y,1);
	  FA[m][k]=*K;
  }
}

  double b=(*Y)(2,0);
  cout << a << endl;
  cout << "hier" << endl;
  double m = FA[1,2]->value(1,1);
  cout <<m << endl;
  FeatureMatrix *FM= new FeatureMatrix(*X,*Y,1);
  Feature **MA;
  MA= FM->FA;
  double n= MA[1,2]->value(1,1);
  MA[1,2]->setlambda(2,2,1.5);
  double t= MA[1,2]->getlambda(2,2);
  cout << n << endl;
  cout << t << endl;
  cout << MA[1,2] << endl;

}
