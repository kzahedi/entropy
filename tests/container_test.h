
#include <cppunit/extensions/HelperMacros.h>

class containerTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(containerTest);
  CPPUNIT_TEST(testFilling);
  CPPUNIT_TEST(testDropping);
  CPPUNIT_TEST(testUniformDiscretisation);
  CPPUNIT_TEST(testUniformDiscretisationUnary);
  CPPUNIT_TEST(testCopy);
  CPPUNIT_TEST(testMax);
  CPPUNIT_TEST(testMin);
  CPPUNIT_TEST(testCopyFunc);
  CPPUNIT_TEST(testExtractColumns);
  CPPUNIT_TEST(testNormaliseColumn);
  CPPUNIT_TEST(testMerge);
  CPPUNIT_TEST(testFillMode);
  CPPUNIT_TEST(testUnique1);
  CPPUNIT_TEST(testUnique1);
  CPPUNIT_TEST(testFind1);
  CPPUNIT_TEST(testFind2);
  CPPUNIT_TEST(testFindList1);
  CPPUNIT_TEST(testFindList2);
  CPPUNIT_TEST(testFind1ByContainer);
  CPPUNIT_TEST(testFind2ByContainer);
  CPPUNIT_TEST(testFindList1ByContainer);
  CPPUNIT_TEST(testFindList2ByContainer);
  CPPUNIT_TEST(testFindList3ByContainer);

  CPPUNIT_TEST_SUITE_END();

  public:

    void testFilling();
    void testDropping();
    void testCopy();
    void testUniformDiscretisation();
    void testUniformDiscretisationUnary();
    void testMax();
    void testMin();
    void testExtractColumns();
    void testNormaliseColumn();
    void testCopyFunc();
    void testMerge();
    void testFillMode();
    void testUnique1();
    void testUnique2();
    void testFind1();
    void testFind2();
    void testFindList1();
    void testFindList2();
    void testFind1ByContainer();
    void testFind2ByContainer();
    void testFindList1ByContainer();
    void testFindList2ByContainer();
    void testFindList3ByContainer();

};
