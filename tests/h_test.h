#include <cppunit/extensions/HelperMacros.h>

class hTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(hTest);
  CPPUNIT_TEST(testMax);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testMax();
    void testZero();

};
