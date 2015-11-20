
#include <cppunit/extensions/HelperMacros.h>

class miTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(miTest);
  CPPUNIT_TEST(testSinus);
  CPPUNIT_TEST(testSparseVsNonSparse);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testSinus();
    void testSparseVsNonSparse();

};
