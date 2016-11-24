#include <cppunit/extensions/HelperMacros.h>

class itMatrixTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(itMatrixTest);
  CPPUNIT_TEST(testFillXBinary);
  CPPUNIT_TEST(testFillXMore);
  CPPUNIT_TEST_SUITE_END();

  public:
    void testFillXBinary();
    void testFillXMore();
};

