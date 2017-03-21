#include <cppunit/extensions/HelperMacros.h>

class deltaTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(deltaTest);
  CPPUNIT_TEST(testDeltaMatch);
  CPPUNIT_TEST(testDeltaMatchXY);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testDeltaMatch();
    void testDeltaMatchXY();

};
