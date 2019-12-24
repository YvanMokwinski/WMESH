#ifndef __header_Exception_h__
#define __header_Exception_h__
#include <exception>

#include <sstream>
#define __FILETRACE__ string("at line ") + dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << __LINE__ ) ).str() + string(" in file ") + string(__FILE__)

#define __TRACE__ string("in function '") + string(__FUNCTION__) + string("' at line ") + dynamic_cast< std::ostringstream & >(( std::ostringstream() << std::dec << __LINE__ ) ).str() + string(" in file ") + string(__FILE__)


class Exception : public exception
{  
 protected:
  string m_message;

 public:
 Exception(string message_) : m_message(message_)
  {
  };

  virtual ~Exception() throw() 
    {
    };
  
  virtual const char* what() const throw()
  {
    return this->m_message.c_str();
  }
  
  string ToString() const 
  {
    return this->m_message;
  };  
};


#endif
