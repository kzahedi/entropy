#include <cppunit/extensions/HelperMacros.h>

class gisTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(gisTest);
  CPPUNIT_TEST(testAND);
  CPPUNIT_TEST(testOR);
  CPPUNIT_TEST(testXOR);
  // CPPUNIT_TEST(testUnique);
  // CPPUNIT_TEST(testUnique2);
  // CPPUNIT_TEST(testMC_W);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testAND();
    void testOR();
    void testXOR();

    void testUnique();
    void testUnique2();
    void testMC_W();

};
