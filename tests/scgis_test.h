#include <cppunit/extensions/HelperMacros.h>

class scgisTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(scgisTest);
  CPPUNIT_TEST(testAND);
  CPPUNIT_TEST(testOR);
  CPPUNIT_TEST(testXOR);
  // CPPUNIT_TEST(testUnique);
  // CPPUNIT_TEST(testUnique2);
  CPPUNIT_TEST(testMC_W);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testAND();
    void testOR();
    void testXOR();

    void testUnique();
    void testUnique2();
    void testMC_W();

};
