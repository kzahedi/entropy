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
void originalTest::testOriginal(){

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
}
