
#include <cppunit/extensions/HelperMacros.h>

class containerTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(containerTest);
  CPPUNIT_TEST(testFilling);
  CPPUNIT_TEST(testDropping);
  CPPUNIT_TEST(testUniformDiscretisation);
  CPPUNIT_TEST(testUniformDiscretisationUnary);
  CPPUNIT_TEST_SUITE_END();

  public:

    void testFilling();
    void testDropping();
    void testUniformDiscretisation();
    void testUniformDiscretisationUnary();

};
