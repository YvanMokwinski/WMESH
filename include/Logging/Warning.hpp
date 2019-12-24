#pragma once
#include "Base.hpp"

namespace Logging
{

#if 0
  class Warning : public Logging::Base<Warning>
  {
    
  public:
    Warning(const int& ithread_)    : Logging::Base<Warning>(ithread_,stdout)
    {

    };
    virtual ~Warning(){};
    static const char *GetToken()
    {
      return "WARNING";
    };
    FILE *GetOutput()const
    {
      return stdout;
    };

  };
#endif
  class Warning : public Logging::Base
  {
    
  public:
    Warning(const int& ithread_)    : Logging::Base(ithread_,stdout,"WARNING")
    {

    };
    virtual ~Warning(){};


    static Warning*Create(const char*		progname_,
			   const unsigned char 	mode_ = STD | LOG)
    {
      const int cpuId 	= Mpi::GetCpuId();
      Warning * m 	= Multiton<Warning,int>::GetInstance(cpuId);
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

    static Warning&Get()
    {
      return *Multiton<Warning,int>::GetInstance(Mpi::GetCpuId());
    };


    FILE *GetOutput()const
    {
      return stdout;
    };

  };
  
};

//#define Warning Warning::Get() 
