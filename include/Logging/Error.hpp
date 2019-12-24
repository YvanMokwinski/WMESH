#pragma once
#include "Base.hpp"

namespace Logging
{
#if 0
  class Error : public Logging::Base<Error>
  {
    
  public:
    Error(const int& ithread_) : Logging::Base<Error>(ithread_,
						      stderr)
    {
    };
    
    virtual ~Error()
    {
      
    };

    static const char *GetToken()
    {
      return "ERROR";
    };
    
  };

  

#endif


    class Error : public Logging::Base
  {
    
  public:
    Error(const int& ithread_) : Logging::Base(ithread_,
					       stderr,"ERROR")
    {
    };
    
    virtual ~Error()
    {
      
    };


    static Error*Create(const char*		progname_,
			   const unsigned char 	mode_ = STD | LOG)
    {
      const int cpuId 	= Mpi::GetCpuId();
      Error * m 	= Multiton<Error,int>::GetInstance(cpuId);
      sprintf	(m->m_progname,"%s",progname_);
      m->m_mode	= mode_;
      unsigned char mode = m->m_mode;
      mode = mode >> 1;
      if (mode%2>0)
	{
	  char a[256];
	  const char * progname_without_path 	= NULL;
	  const char * tmp 			= NULL;
	  for (tmp=m->m_progname;tmp[0]!='\0';++tmp)
	    if (*tmp=='/')
	      progname_without_path = tmp;
	  progname_without_path = (progname_without_path)?progname_without_path+1:m->m_progname;
	  sprintf(a,"%s.cpu%d.log",
		  progname_without_path,
		  cpuId);	  
	  m->m_logfile = fopen(a,"w");
	}
      return m;
    };

    static Error&Get()
    {
      return *Multiton<Error,int>::GetInstance(Mpi::GetCpuId());
    };

    
    static const char *GetToken()
    {
      return "ERROR";
    };
    
  };

  

};

//#define Error Error::Get() 
