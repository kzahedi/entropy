#include <cppunit/extensions/HelperMacros.h>

class originalTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(originalTest);
   CPPUNIT_TEST(testOriginal);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testOriginal();
};
