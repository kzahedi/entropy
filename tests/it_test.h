#include <cppunit/extensions/HelperMacros.h>

class itTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(itTest);
  CPPUNIT_TEST(OneXOneY);
  CPPUNIT_TEST(TwoXOneY);
  CPPUNIT_TEST(TwoXTwoY);
  CPPUNIT_TEST_SUITE_END();

  public:

    void OneXOneY();
    void TwoXOneY();
    void TwoXTwoY();

};
