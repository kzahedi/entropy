
#include <cppunit/extensions/HelperMacros.h>

class cmiTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(cmiTest);
  CPPUNIT_TEST(testSinus);
  CPPUNIT_TEST(testSparseVsNonSparse);
  CPPUNIT_TEST(testMatrixWiseComparision);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testSinus();
    void testSparseVsNonSparse();
    void testMatrixWiseComparision();

};
