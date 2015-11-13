#ifndef __MATRIX_EXCEPTION_H__
#define __MATRIX_EXCEPTION_H__

#include <vector>
#include <string>
#include <iostream>

class MatrixException : public std::exception
{
  public:
    explicit MatrixException(const std::string& what)
      :m_what(what)
    {}

    virtual ~MatrixException() throw() {}

    virtual const char * what() const throw()
    {
      return m_what.c_str();
    }

    virtual void message() const throw()
    {
      std::cout << "MatrixException: " << m_what << std::endl;
    }


  private:
    std::string m_what;
};

#endif // __MATRIX_EXCEPTION_H__
