#include <cppunit/extensions/HelperMacros.h>

class iterativeScalingTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(iterativeScalingTest);
  CPPUNIT_TEST(ValueTest);
  // CPPUNIT_TEST(OneXOneY);
  // CPPUNIT_TEST(OneXOneYFour);
  // CPPUNIT_TEST(OneXOneYEight);
  // CPPUNIT_TEST(OneXOneYTwenty);
  // CPPUNIT_TEST(TwoXOneY);
  // CPPUNIT_TEST(TwoXTwoY);
  // CPPUNIT_TEST(NotBinarySix);
  CPPUNIT_TEST_SUITE_END();


  public:
    void ValueTest();
    // void OneXOneY();
    // void OneXOneYFour();
    // void OneXOneYEight();
    // void OneXOneYTwenty();
    // void TwoXOneY();
    // void TwoXTwoY();
    // void NotBinarySix();
};

