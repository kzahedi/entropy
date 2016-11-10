#include <cppunit/extensions/HelperMacros.h>

class itTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(itTest);
  CPPUNIT_TEST(ValueTest);
  CPPUNIT_TEST(OneXOneY);
  CPPUNIT_TEST(OneXOneYFour);
  CPPUNIT_TEST(OneXOneYEight);
  CPPUNIT_TEST(OneXOneYTwenty);
  CPPUNIT_TEST(TwoXOneY);
  CPPUNIT_TEST(TwoXTwoY);
  CPPUNIT_TEST(NotBinarySix);
//  CPPUNIT_TEST(FourXFourY);
  CPPUNIT_TEST_SUITE_END();


  public:
    void ValueTest();
    void OneXOneY();
    void OneXOneYFour();
    void OneXOneYEight();
    void OneXOneYTwenty();
    void TwoXOneY();
    void TwoXTwoY();
    void NotBinarySix();
  //  void FourXFourY();
};

