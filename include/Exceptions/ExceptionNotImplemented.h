#ifndef __header_ExceptionNotImplemented_h__
#define __header_ExceptionNotImplemented_h__

#include "Exception.h"

class ExceptionNotImplemented : public Exception
{  
public:
 ExceptionNotImplemented(string message_) : Exception(message_)
  {
  };
};


#endif
