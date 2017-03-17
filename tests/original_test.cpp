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
  cout << endl;
  double m0_0 = Test->getMarginalProp(0,2,p);
  double m0_1 = Test->getMarginalProp(7,2,p);
  double m1_0 = Test->getMarginalProp(0,1,p);
  double m1_1 = Test->getMarginalProp(7,1,p);
  double m2_0 = Test->getMarginalProp(0,0,p);
  double m2_1 = Test->getMarginalProp(7,0,p);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[0], m0_0*m1_0*m2_0,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[1], m0_0*m1_0*m2_1,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[2], m0_0*m1_1*m2_0,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[3], m0_0*m1_1*m2_1,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[4], m0_1*m1_0*m2_0,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[5], m0_1*m1_0*m2_1,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[6], m0_1*m1_1*m2_0,0.0001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(pconv[7], m0_1*m1_1*m2_1,0.0001);
//  cout << m0_0 << " "<<  m0_1 << " " << m1_0 << " " << m1_1 << " "<< m2_0 << " " << m2_1 << " " << endl;
}
