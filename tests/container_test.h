#ifndef __CONFIGURATION_TEST_H__
#define __CONFIGURATION_TEST_H__

#include <cppunit/extensions/HelperMacros.h>

class containerTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(containerTest);
  CPPUNIT_TEST(testFilling);
  CPPUNIT_TEST(testUniformDiscretisation);
  CPPUNIT_TEST(testUniformDiscretisationUnary);
  CPPUNIT_TEST_SUITE_END();

  public:

  void testFilling();
  void testUniformDiscretisation();
  void testUniformDiscretisationUnary();

};

#endif // __RNN_H__
