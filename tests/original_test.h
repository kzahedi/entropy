#include <cppunit/extensions/HelperMacros.h>

class originalTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(originalTest);
   CPPUNIT_TEST(marginalFeatures);
   CPPUNIT_TEST(neighbourhoodRelations);
   CPPUNIT_TEST(testXor);
   CPPUNIT_TEST(testOr);
   CPPUNIT_TEST(testAnd);
  CPPUNIT_TEST_SUITE_END();

  public:

    void marginalFeatures();
    void neighbourhoodRelations();
    void testXor();
    void testOr();
    void testAnd();
};
