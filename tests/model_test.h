#include <cppunit/extensions/HelperMacros.h>

class modelTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(modelTest);
  CPPUNIT_TEST(testUnique);
  // CPPUNIT_TEST(testUnique2);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testUnique();
    void testUnique2();

};
