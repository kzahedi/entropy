#include <cppunit/extensions/HelperMacros.h>

class entropyTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(entropyTest);
  CPPUNIT_TEST(testMax);
  CPPUNIT_TEST(testZero);
  CPPUNIT_TEST(testConditional);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testMax();
    void testZero();
    void testConditional();

};
