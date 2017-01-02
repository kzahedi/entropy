#include <cppunit/extensions/HelperMacros.h>

class isoTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(isoTest);
  CPPUNIT_TEST(test);
  CPPUNIT_TEST_SUITE_END();

  public:

    void test();

};
