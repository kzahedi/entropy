#include "iterative_scaling_test.h"

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/sparse/MI.h>

#include <entropy++/iterativescaling/IterativeScalingBase.h>
#include <entropy++/iterativescaling/gis/IterativeScaling.h>
#include <entropy++/iterativescaling/scgis/IterativeScaling.h>

#include <iostream>
#include <string>

using namespace std;
using namespace entropy;
using namespace entropy::iterativescaling;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(iterativeScalingTest);

void iterativeScalingTest::ValueTest()
{
#ifdef __APPLE__
  sranddev();
#else
  srand(time(NULL));
#endif

  ULContainer *xAlphabet = new ULContainer(2,1); // alphabet
  *xAlphabet << 0 << 1;
  ULContainer *yAlphabet = new ULContainer(2,1); // alphabet
  *yAlphabet << 0 << 1;

  ivvector alphX(1,ivector(0));
  alphX[0].push_back(0);

  ivvector alphY(1,ivector(0));
  alphY[0].push_back(0);

  dvector lambda(1);
  lambda[0] = 1;

  int xColumns      = 2;
  int xRows         = 1000;
  int yColumns      = 1;

  ULContainer* xData = new ULContainer(xRows,xColumns);

  for(int i = 0;i < xRows; i++)
  {
    for(int j = 0;j < xColumns; j++)
    {
      double z = rand() % xAlphabet->rows(); // random indices for the alphabet
      *xData << (*xAlphabet)(z,0);
    }
  }

  IterativeScalingBase *exact = new IterativeScalingBase(yColumns,
                                                         xData,
                                                         xAlphabet,
                                                         yAlphabet,
                                                         alphX, // TODO what is alphX/Y
                                                         alphY);

  // for(int feat = 0; feat < alphX.size(); feat++)
  // {
    // for(int delti = 0; delti < pow(xAlphabet->rows(), alphX[feat].size()); delti++)
    // {
      // for(int deltj = 0; deltj < pow(yAlphabet->rows(), alphY[feat].size()); deltj++)
      // {
        // double r = rand() % lambda.size();
        // exact->setFeatureArraylambda(feat, delti, deltj, lambda[r]);
      // }
    // }
  // }

  // IContainer* generatedY = new IContainer(xRows, yColumns);
  // for(int i = 0; i < xRows; i++)
  // {
    // double z   = (double)rand()/RAND_MAX;
    // double s   = 0;
    // int    ind = 0;
    // for(int j = 0; j < pow(yAlphabet->rows(), yColumns) && s < z; j++)
    // {
      // s   += exact->prop(i,j);
      // ind  = j;
    // }
    // for(int k = 0; k < yColumns; k++)
    // {
      // (*generatedY) << exact->index(ind, false, yColumns)[k];
    // }
  // }


  // Test *test_1 = new Test(1,1,100,lambda, *zX,*zY,alphX, alphY);
  for(int x = 0; x < 2; x++)
  {
    for(int y = 0; y < 2; y++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, exact->p_x_c_y(x, y), 0.1);
    }
  }

  IContainer *indizes = new IContainer(1,3);
  (*indizes) << 0 << 0 << 0;
  DContainer *values  = new DContainer(1,1);
  (*values) << 5;

  for(int i = 0; i < indizes->rows(); i++)
  {
    exact->setFeatureArraylambda(indizes->get(i,0), indizes->get(i,1), indizes->get(i,2), values->get(i,0));
  }




  CPPUNIT_ASSERT(exact->p_x_c_y(0,0) > exact->p_x_c_y(0,1));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, exact->p_x_c_y(1,1), 0.01);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, exact->p_x_c_y(1,0), 0.01);

}

// void iterativeScalingTest::OneXOneY()
// {
  // cout << "OneXOneY" << endl;
  // srand(time(NULL));
  // DContainer *zX = new DContainer(2,1); // alphabet
  // *zX << 0 << 1;
  // DContainer *zY = new DContainer(2,1); // alphabet
  // *zY << 0 << 1;
  // dvector lambda(3);
  // lambda[0] = 0;
  // lambda[1] = 1;
  // lambda[2] = 5;

  // vector<vector<int > > alphX(1,ivector(0));
  // alphX[0].push_back(0);

  // vector<vector<int > > alphY(1,ivector(0));
  // alphY[0].push_back(0);

  // IsParameter param;
  // param.lambdavalue    = 1.0;
  // param.lambdadeltaval = 1.0;
  // param.sigma          = 0.01;
  // param.maxit          = 500;
  // param.konv           = 0.00001;
  // param.konvtime       = true;
  // param.time           = false;
  // param.test           = true;
  // param.seconds        = 72000;

  // Test *testval = new Test(1,1,10000,lambda, *zX,*zY,alphX, alphY); //get data

  // GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zX,*zY,alphX,alphY,param);
  // SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zX,*zY,alphX,alphY,param);

  // p(y_0 = 0 | x_0 = 0) + p(y_0 = 1 | x_0 = 0) = 1.0
  // CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->propAlphX( 0, 0)+test->propAlphX( 0, 1), 0.1);
  // p(y_0 = 0 | x_0 = 1) + p(y_0 = 1 | x_0 = 1) = 1.0
  // CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->propAlphX( 1, 0)+test->propAlphX( 1, 1), 0.1);

  // CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->propAlphX(0, 0)+testsc->propAlphX(0, 1), 0.1);
  // CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->propAlphX(1, 0)+testsc->propAlphX(1, 1), 0.1);

  // double exact;
  // for(int i=0;i< 2; i++)
  // {
    // for(int j=0; j< 2; j++)
    // {
      // exact= testval->getProp(i,j);
      // cout << i << " " << j << endl;
      // cout << "exakt: " <<  testval->getProp(i,j) << endl;
      // cout << "gis:  " << test->propAlphX(i,j) << endl;
      // cout << "scgis: " << testsc->propAlphX(i,j) << endl;
    // }
  // }

  // getsizeconv -> size of stored l-values (constraints)
  // constraints: difference betweens observed feature frequency vs. expected feature frequency, given the lambdas
  // for(int i = 0; i < test->getsizeconv()-1; i++)
  // {
    // if((test->getconv(i)>pow(10,-11)))
    // {
      // if(!(test->getconv(i) >= test->getconv(i+1))) // must be decreasing, for output
      // {
        // cout << test->getconv(i) << " " << test->getconv(i+1) << endl;
      // }
      // CPPUNIT_ASSERT(test->getconv(i) >= test->getconv(i+1));
    // }
  // }
  // for(int i=0; i< testsc->getsizeconv()-1; i++ )
  // {
    // if((testsc->getconv(i)>pow(10,-11)))
    // {
      // if(!(testsc->getconv(i) >= testsc->getconv(i+1)))
      // {
        // cout << testsc->getconv(i)  << " " << testsc->getconv(i+1) << endl;
      // }
      // CPPUNIT_ASSERT  (testsc->getconv(i) >= testsc->getconv(i+1) );
    // }
  // }
  // cout << " GIS: " << testval->KL(test);
  // cout << " SCGIS: " << testval->KL(testsc);
// }

