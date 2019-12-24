#ifndef __header_ExceptionDebug_h__
#define __header_ExceptionDebug_h__
#include "Exception.h"

class ExceptionDebug : public Exception
{  
public:
 ExceptionDebug(string message_) : Exception(message_)
  {
  };
};


#endif
