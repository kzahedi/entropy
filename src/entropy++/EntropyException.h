#ifndef __ENTROPY_EXCEPTION_H__
#define __ENTROPY_EXCEPTION_H__

#include <vector>
#include <string>
#include <iostream>

class EntropyException : public std::exception
{
  public:
    explicit EntropyException(const std::string& what)
      :m_what(what)
    {}

    virtual ~EntropyException() throw() {}

    virtual const char * what() const throw()
    {
      return m_what.c_str();
    }

    virtual void message() const throw()
    {
      std::cout << "EntropyException: " << m_what << std::endl;
    }


  private:
    std::string m_what;
};

#endif // __ENTROPY_EXCEPTION_H__
