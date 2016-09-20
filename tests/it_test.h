#include <cppunit/extensions/HelperMacros.h>

class itTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(itTest);
  CPPUNIT_TEST(OneXOneY);
  CPPUNIT_TEST(SCOneXOneY);
  CPPUNIT_TEST(TwoXOneY);
  CPPUNIT_TEST(SCTwoXOneY);
  CPPUNIT_TEST(TwoXTwoY);
  CPPUNIT_TEST(SCTwoXTwoY);
  CPPUNIT_TEST(NotBinary);
  CPPUNIT_TEST(FourXFourY);
  CPPUNIT_TEST_SUITE_END();


  public:
    void OneXOneY();
    void SCOneXOneY();
    void TwoXOneY();
    void SCTwoXOneY();
    void TwoXTwoY();
    void SCTwoXTwoY();
    void NotBinary();
    void FourXFourY();

};

