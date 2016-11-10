#include <cppunit/extensions/HelperMacros.h>

class csvTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(csvTest);
  /*
  CPPUNIT_TEST(readTestFile);
  */
  CPPUNIT_TEST_SUITE_END();

  public:

    void readTestFile();

};
