#include <cppunit/extensions/HelperMacros.h>

class matrixTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(matrixTest);
  CPPUNIT_TEST(testInitialisation);
  CPPUNIT_TEST(testSet);
  CPPUNIT_TEST(testAdd);
  CPPUNIT_TEST(testMul);
  CPPUNIT_TEST(testDiv);
  CPPUNIT_TEST_SUITE_END();

  public:

  void testInitialisation();
  void testSet();
  void testAdd();
  void testMul();
  void testDiv();

};




// SparseMatrix&   operator-=(const SparseMatrix &m)  throw(MatrixException);

// const SparseMatrix operator* (const double factor);
// SparseMatrix&      operator*=(const double factor);

// SparseMatrix&      operator= (const SparseMatrix &m);

// SparseMatrix&      operator-=(const double d);
// const SparseMatrix operator-(const SparseMatrix &m);

