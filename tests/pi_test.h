
#include <cppunit/extensions/HelperMacros.h>

class piTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(piTest);
  CPPUNIT_TEST(testSinus);
  CPPUNIT_TEST(testSparseVsNonSparse);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testSinus();
    void testSparseVsNonSparse();

};
