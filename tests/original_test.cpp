#include "original_test.h"

#include <entropy++/Matrix.h>
#include <entropy++/iterativescaling/Original.h>

#include <iostream>
#include <string>

using namespace std;
using namespace entropy;
using namespace entropy::iterativescaling;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( originalTest );
void originalTest::marginalFeatures(){

  vector<vector<int> > features;

  vector<int> a;
  a.push_back(0);
  features.push_back(a);
  vector<int> b;
  b.push_back(1);
  features.push_back(b);
  vector<int> c;
  c.push_back(2);
  features.push_back(c);

  vector<double> p = vector<double> (8);
  p[0]=0.5;
  p[4]=0.25;
  p[6]=0.25;

  Original* Test = new Original(3, features, p );
  Test->iterate(0.001);
  vector<double> pconv = Test->getp();
//  for(int i=0;i<8;i++){
//    cout<< pconv[i] <<"  ";
//  }
///  cout << endl;
  vector<int> m1;
  m1.push_back(0);
  vector<int> m2;
  m2.push_back(1);
  vector<int> m3;
  m3.push_back(2);
  double m0_0 = Test->getMarginalProp(0,m1,p);
  double m0_1 = Test->getMarginalProp(7,m1,p);
  double m1_0 = Test->getMarginalProp(0,m2,p);
  double m1_1 = Test->getMarginalProp(7,m2,p);
  double m2_0 = Test->getMarginalProp(0,m3,p);
  double m2_1 = Test->getMarginalProp(7,m3,p);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[0], m0_0*m1_0*m2_0,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[1], m0_0*m1_0*m2_1,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[2], m0_0*m1_1*m2_0,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[3], m0_0*m1_1*m2_1,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[4], m0_1*m1_0*m2_0,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[5], m0_1*m1_0*m2_1,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[6], m0_1*m1_1*m2_0,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[7], m0_1*m1_1*m2_1,0.0001);
// cout << m0_0 << " "<<  m0_1 << " " << m1_0 << " " << m1_1 << " "<< m2_0 << " " << m2_1 << " " << endl;
//  cout << (*Test) << endl;
}
void originalTest:: neighbourhoodRelations(){
  vector<vector<int> > features;

  vector<int> a;
  a.push_back(0);
  a.push_back(1);
  features.push_back(a);

  vector<int> b;
  b.push_back(1);
  b.push_back(2);
  features.push_back(b);

  vector<int> c;
  c.push_back(2);
  c.push_back(3);
  features.push_back(c);

  vector<double> p = vector<double> (16);
    p[0]=0.15;
    p[1]=0.05;
    p[4]=0.25;
    p[6]=0.15;
    p[9]=0.1;
    p[12]=0.3;

  Original* Test = new Original(4, features, p );
  Test->iterate(0.001);
  vector<double> pconv = Test->getp();
 // for(int i=0;i<8;i++){
 //   cout<< pconv[i] << endl;
 //  }
  vector<int> m1;
  m1.push_back(0);
  vector<int> m2;
  m2.push_back(1);
  vector<int> m3;
  m3.push_back(2);
  vector<int> m4;
  m4.push_back(3);

  double prop =0.0;
  for(int i=0;i<16;i++){
    prop = Test->getMarginalProp(i,m1,p)*Test->getConditionalProp(m2,m1,i,p)*Test->getConditionalProp(m3,m2,i,p)*Test->getConditionalProp(m4,m3,i,p);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(prop, pconv[i],0.0001);
   // cout <<Test->getMarginalProp(i,m1,p)<< " " << Test->getConditionalProp(m2,m1,i,p) << " " << Test->getConditionalProp(m3,m2,i,p)<< " "<< Test->getConditionalProp(m4,m3,i,p)<<" "<<   prop<<" "<<pconv[i] << endl;
  }
 // cout << (*Test) << endl;

}
void originalTest:: testXor(){
  vector<vector<int> > features;

  vector<int> a;
  a.push_back(0);
  a.push_back(2);
  features.push_back(a);

  vector<int> b;
  b.push_back(1);
  b.push_back(2);
  features.push_back(b);

  vector<double> p = vector<double> (8);
  p[0]=0.25;
  p[3]=0.25;
  p[5]=0.25;
  p[6]=0.25;

  Original* Test = new Original(3, features, p );
   Test->iterate(5);
   vector<double> pconv = Test->getp();
 //  for(int i=0;i<8;i++){
 //    cout<< pconv[i] << endl;
//   }
   vector<int> x;
    x.push_back(0);
    x.push_back(1);
  vector<int> y;
    y.push_back(2);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(Test->calculateKL(p,pconv), 1 ,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(Test->calculateConditionalKL(p,pconv,y,x) , 1 ,0.0001);
  //cout << Test->calculateKL(p,pconv)<< endl;
  //cout << Test->calculateConditionalKL(p,pconv,y,x)  << endl;
}
void originalTest:: testOr(){
  vector<vector<int> > features;

  vector<int> a;
  a.push_back(0);
  a.push_back(2);
  features.push_back(a);

  vector<int> b;
  b.push_back(1);
  b.push_back(2);
  features.push_back(b);

  vector<double> p = vector<double> (8);
  p[0]=0.25;
  p[3]=0.25;
  p[5]=0.25;
  p[7]=0.25;

  Original* Test = new Original(3, features, p );
   Test->iterate(0.0001);
   vector<double> pconv = Test->getp();
 //  for(int i=0;i<8;i++){
 //    cout<< pconv[i] << endl;
//   }
   vector<int> x;
    x.push_back(0);
    x.push_back(1);
  vector<int> y;
    y.push_back(2);

//  CPPUNIT_ASSERT_DOUBLES_EQUAL(Test->calculateKL(p,pconv), 0.5,0.0001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Test->calculateKL(p,pconv), 0.188722 ,0.0001);
//  cout << Test->calculateKL(p,pconv)<< endl;
 //  cout << Test->calculateConditionalKL(p,pconv,y,x)  << endl;
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Test->calculateConditionalKL(p,pconv,y,x), 0.103759, 0.0001);
}
void originalTest:: testAnd(){
  vector<vector<int> > features;

  vector<int> a;
  a.push_back(0);
  a.push_back(2);
  features.push_back(a);

  vector<int> b;
  b.push_back(1);
  b.push_back(2);
  features.push_back(b);

  vector<double> p = vector<double> (8);
  p[0]=0.25;
  p[2]=0.25;
  p[4]=0.25;
  p[7]=0.25;

  Original* Test = new Original(3, features, p );
   Test->iterate(0.0001);
   vector<double> pconv = Test->getp();
 //  for(int i=0;i<8;i++){
 //    cout<< pconv[i] << endl;
//   }
   vector<int> x;
    x.push_back(0);
    x.push_back(1);
   vector<int> y;
    y.push_back(2);


//  CPPUNIT_ASSERT_DOUBLES_EQUAL(Test->calculateKL(p,pconv), 0.5,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(Test->calculateKL(p,pconv), 0.188722 ,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(Test->calculateConditionalKL(p,pconv,y,x), 0.103759, 0.0001);
 // cout << Test->calculateKL(p,pconv)<< endl;
 // cout << Test->calculateConditionalKL(p,pconv,y,x)  << endl;
}
