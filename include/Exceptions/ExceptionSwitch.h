#ifndef __header_ExceptionSwitch_h__
#define __header_ExceptionSwitch_h__
#include "Exception.h"

class ExceptionSwitch : public Exception
{  
public:
 ExceptionSwitch(string message_) : Exception(message_)
    {
    };
};

#define SWITCH_EXCEPTION throw(ExceptionSwitch("switch failed"))


#endif
