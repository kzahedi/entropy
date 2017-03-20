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
  cout << "vor Original" << endl;
  Original* Test = new Original(3, features, p );
  Test->iterate(6);
  vector<double> pconv = Test->getp();
  for(int i=0;i<8;i++){
    cout<< pconv[i] <<"  ";
  }
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

  vector<double> p = vector<double> (8);
    p[0]=0.45;
    p[1]=0.05;
    p[4]=0.25;
    p[6]=0.25;

  Original* Test = new Original(3, features, p );
  Test->iterate(10);
  vector<double> pconv = Test->getp();
  for(int i=0;i<8;i++){
    cout<< pconv[i] << endl;
  }
  vector<int> m1;
  m1.push_back(0);
  vector<int> m2;
  m2.push_back(1);
  vector<int> m3;
  m3.push_back(2);
  double prop =0.0;
  for(int i=0;i<8;i++){
    prop = Test->getMarginalProp(i,m1,p)*Test->getConditionalProp(m2,m1,i,p)*Test->getConditionalProp(m3,m2,i,p);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(prop, pconv[i],0.0001);
    cout <<Test->getMarginalProp(i,m1,p)<< " " << Test->getConditionalProp(m2,m1,i,p) << " " << Test->getConditionalProp(m3,m2,i,p)<< " " << prop << endl;
  }

}
