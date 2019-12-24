#ifndef __header_ExceptionAssert_h__
#define __header_ExceptionAssert_h__
#include "Exception.h"

class ExceptionAssert : public Exception
{  
public:
 ExceptionAssert(string message_) : Exception(message_)
  {
  };
};


#endif
