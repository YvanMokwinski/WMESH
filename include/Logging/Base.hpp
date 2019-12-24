#pragma once

//#include "System/System.h"
// #include "System.h"
#include "Multiton.h"
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <valarray>


namespace Logging
{  

  class Base : public Multiton<Base,int>
  {
  public:
    /**
       \brief Define mode for monitoring   
       @note an example should be a=LoggingMode_STD | LoggingMode_LOG
    */  
    typedef enum __eLoggingMode
      { 
	OFF		= 0,/*!< default logical value*/
	STD		= 1,/*!< mode to display messages on standard output*/
	LOG		= 2,/*!< mode to display messages on log output*/
	EXITIFERROR	= 4 /*!< mode to exit if an error message is received */
      } Mode;

  private:
    friend Base* 	Multiton<Base,int>::GetInstance(const int&);
    friend void 	Multiton<Base,int>::KillInstances();
  protected:    
    //  string 		m_disp;
    std::stringstream 	m_s;
    FILE * 		m_output;
    unsigned char 	m_mode;
    char 		m_progname[256];
    FILE * 		m_logfile;
    int 		m_ithread;
    const char * m_token;
  protected:
    
    Base(const int& ithread_,	 
	 FILE * file_,
	 const char * token_) 
    {
      m_token = token_;
      m_output = file_;
      //      m_disp = disp_;
      char	ctmp[256];
      m_ithread = ithread_;
      int len = Mpi::GetCpuName(m_ithread,ctmp);
      len += sprintf(&ctmp[len],".%s:",token_);   
      m_s << ctmp;
    };
  
    
 public:
    
    template<typename T> Base& operator<<(const T& v) 
    {
      m_s << v;
      return *this;      
    };
    
    Base const& operator<<(std::ostream& (*F)(std::ostream&))
    { 
      char	ctmp[256];
      unsigned char mode = m_mode;   
      if (mode%2>0)
	{
	  /*LoggingMode_STD*/
	  fflush(m_output);
	  fprintf(m_output,"//%s\n",m_s.str().c_str());
	}
      mode = mode >> 1;
      if (mode%2>0)
	{
	  /*LoggingMode_LOG*/
	  fprintf(m_logfile,"//%s\n",m_s.str().c_str());
	}
      m_s.str("");
      m_s.clear();
      int len = Mpi::GetCpuName(m_ithread,ctmp);
      len += sprintf(&ctmp[len],".%s:",m_token);   
      m_s << ctmp;
      return *this; 
    };
    

#if 0

        static Base*Create(const char*		progname_,
			   const unsigned char 	mode_ = STD | LOG)
    {
      const int cpuId 	= Mpi::GetCpuId();
      Base * m 	= Multiton<Base,int>::GetInstance(cpuId);
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
    static Base&Get()
    {
      const int cpuId = Mpi::GetCpuId();  
      return *Multiton<Base,int>::GetInstance(cpuId);
    };
#endif
    
    
    /**
       \brief Shutdown the monitor
       @note it will close all log files.
    */
    virtual ~Base()
    {    
      if (m_logfile)
	{
	  std::cerr << "close " << std::endl;
	  fclose(m_logfile);
	}
    };
    
  };

#if 0  
  template <typename tLogging> class Base : public Multiton<Base<tLogging>,int>
  {
  public:
    /**
       \brief Define mode for monitoring   
       @note an example should be a=LoggingMode_STD | LoggingMode_LOG
    */  
    typedef enum __eLoggingMode
      { 
	OFF		= 0,/*!< default logical value*/
	STD		= 1,/*!< mode to display messages on standard output*/
	LOG		= 2,/*!< mode to display messages on log output*/
	EXITIFERROR	= 4 /*!< mode to exit if an error message is received */
      } Mode;

  private:
    friend tLogging* 	Multiton<tLogging,int>::GetInstance(const int&);
    friend void 	Multiton<tLogging,int>::KillInstances();
  protected:    
    //  string 		m_disp;
    std::stringstream 	m_s;
    FILE * 		m_output;
    unsigned char 	m_mode;
    char 		m_progname[256];
    FILE * 		m_logfile;
    int 		m_ithread;
  protected:
    
    Base(const int& ithread_,
	 FILE * file_) 
    {
      m_output = file_;
      //      m_disp = disp_;
      char	ctmp[256];
      m_ithread = ithread_;
      int len = Mpi::GetCpuName(m_ithread,ctmp);
      len += sprintf(&ctmp[len],".%s:",tLogging::GetToken());   
      m_s << ctmp;
    };
  
    
 public:
    
    template<typename T> Base<tLogging>& operator<<(const T& v) 
    {
      m_s << v;
      return *this;      
    };
    
    Base<tLogging> const& operator<<(std::ostream& (*F)(std::ostream&))
    { 
      char	ctmp[256];
      unsigned char mode = m_mode;   
      if (mode%2>0)
	{
	  /*LoggingMode_STD*/
	  fflush(m_output);
	  fprintf(m_output,"//%s\n",m_s.str().c_str());
	}
      mode = mode >> 1;
      if (mode%2>0)
	{
	  /*LoggingMode_LOG*/
	  fprintf(m_logfile,"//%s\n",m_s.str().c_str());
	}
      m_s.str("");
      m_s.clear();
      int len = Mpi::GetCpuName(m_ithread,ctmp);
      len += sprintf(&ctmp[len],".%s:",tLogging::GetToken());   
      m_s << ctmp;
      return *this; 
    };
    

    static tLogging*Create(const char*		progname_,
			   const unsigned char 	mode_ = STD | LOG)
    {
      const int cpuId 	= Mpi::GetCpuId();
      tLogging * m 	= Multiton<tLogging,int>::GetInstance(cpuId);
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
    
    static tLogging&Get()
    {
      const int cpuId = Mpi::GetCpuId();  
      return *Multiton<tLogging,int>::GetInstance(cpuId);
    };
    
    /**
       \brief Shutdown the monitor
       @note it will close all log files.
    */
    virtual ~Base()
    {    
      if (m_logfile)
	{
	  std::cerr << "close " << std::endl;
	  fclose(m_logfile);
	}
    };
    
  };
  
#endif
  class Message : public Base
  {
  public:
    Message(const int& ithread_) : Base(ithread_,stdout,"MSG")
    {
      
    };
    virtual ~Message(){};

    static Message&Get()
    {
      const int cpuId = Mpi::GetCpuId();  
      return *Multiton<Message,int>::GetInstance(cpuId);
    };

    static Message*Create(const char*		progname_,
			   const unsigned char 	mode_ = STD | LOG)
    {
      const int cpuId 	= Mpi::GetCpuId();      
      Message& m 	= Message::Get();
      sprintf	(m.m_progname,"%s",progname_);
      m.m_mode	= mode_;
      unsigned char mode = m.m_mode;
      mode = mode >> 1;
      if (mode%2>0)
	{
	  char a[256];
	  const char * progname_without_path 	= NULL;
	  const char * tmp 			= NULL;
	  for (tmp=m.m_progname;tmp[0]!='\0';++tmp)
	    if (*tmp=='/')
	      progname_without_path = tmp;
	  progname_without_path = (progname_without_path)?progname_without_path+1:m.m_progname;
	  sprintf(a,"%s.cpu%d.log",
		  progname_without_path,
		  cpuId);	  
	  m.m_logfile = fopen(a,"w");
	}
      return &m;
    };

  };

  
#define LogMessage Logging::Message::Get()
  
};


