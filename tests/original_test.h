#include <cppunit/extensions/HelperMacros.h>

class originalTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(originalTest);
   CPPUNIT_TEST(marginalFeatures);
  CPPUNIT_TEST_SUITE_END();

  public:

    void marginalFeatures();
};
