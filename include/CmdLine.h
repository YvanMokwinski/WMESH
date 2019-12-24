#ifndef __header_CmdLine_h__
#define __header_CmdLine_h__
#include <sys/types.h>
#include <unistd.h>

//#include "System/System.h"
//#include "System.h"
#include <vector>
#include <string>
//#include "Output.h"
//!@brief Class to define the command line arguments.
class CmdLine
{
  friend std::ostream& operator<<(std::ostream& out_, const CmdLine& cmdLine_);
 protected:
  int 			argc;
  std::vector<std::string> 	argv;
  int 			ref_argc;
  char**  		ref_argv;  
  
  int search(int 		argc_,
	     const std::string& 	opt_)
  {
    int i = (int)0;
    for (i=1;i<argc_;i++)
      if (this->argv[i]==opt_)
	break;
    return (i<argc_) ?  i :  -1;
  };
  
 public:
  //!@brief Constructor.
  //!@param argc_ number of arguments.
  //!@param argv_ array of arguments.
  CmdLine(const int		argc_,
	  char** 		argv_)
    {  
      this->ref_argc 	= argc_;
      this->ref_argv 	= argv_;
      this->argc	= argc_;
      for (int i=0;i<argc_;++i) 
	this->argv.push_back(std::string(argv_[i]));      
    };
  
  ~CmdLine()
    {      
    };

  
  bool 		IsEmpty	() 	const { return (this->argc==1); };  
  int 		GetNumArgs() 	const { return this->argc; };
    
  std::string 	GetArg(const int i_)
  {
#if 0
#ifndef NDEBUG
    Debug::IsInRange(__TRACE__,i_,0,this->argc-1);
#endif
#endif
    return this->argv[i_];
  };

  
  bool  Get(const std::string& opt_)
  {  
    const int k = this->search(this->argc,opt_);
    if (k==-1)
      return false;
    if (this->argc<=k)
      return false;
    for (int i=0;i<this->argc-k-1;++i)
      {		
	this->argv[k+i] = this->argv[k+1+i];
      }
    this->argc-=1;
    return true;
  };

  bool  Get(const std::string& opt_,bool*value_)
  { 
    *value_ = Get(opt_);
    return value_;
  };
  
  bool  Get(const std::string& opt_,std::string& t_)
  {  
    const int k = this->search(this->argc,opt_);
    if (k==-1)
      return false;
    if (this->argc <= k + 1)
      return false;    
    t_ = this->argv[k+1];
    for (int i=0;i<this->argc-k-2;++i)
      this->argv[k+i] = this->argv[2+k+i];
    this->argc-=2;
    return true;
  };

  bool  Get(const std::string& opt_,std::string* t_)
  {  
    const int k = this->search(this->argc,opt_);
    if (k==-1)
      return false;
    if (this->argc <= k + 1)
      return false;    
    *t_ = this->argv[k+1];
    for (int i=0;i<this->argc-k-2;++i)
      this->argv[k+i] = this->argv[2+k+i];
    this->argc-=2;
    return true;
  };

  template <typename _int_t> bool   Get(const std::string&	opt_,
	     _int_t* 		x_,
	     const int nb_=1)
  {  
    int k = this->search(this->argc,opt_);
    if (k==-1)
      return false;
    if (this->argc <= k + 1)
      return false;
    for (int i=0;i<nb_;++i)
      sscanf(this->argv[k+1+i].c_str(),"%lld",&x_[i]);
    
    for ( int i=0;i<this->argc-k-1-nb_;++i)
      this->argv[k+i] = this->argv[1+nb_+k+i];
    this->argc-=1+nb_;
    return true;
    
  };


  bool   Get(const std::string&	opt_,
	     double* 			x_)
  {  
    int k = this->search(this->argc,
			 opt_);
    if (k==-1)
      return false;
    if (this->argc <= k + 1)
      return false;
    sscanf(this->argv[k+1].c_str(),"%le",x_);
    for ( int i=0;i<this->argc-k-2;++i)
      this->argv[k+i] = this->argv[2+k+i];
    this->argc-=2;
    return true;
    
  };
  
  void CheckInvalid()const 
  {
#if 0
    for (int i=1;i<this->argc;++i)
      if (this->argv[i][0]=='-')
	MNS_THROW(std::string("Invalid argument '") + this->argv[i] + std::string("'"),__eErr_user);
#endif
  };
 


};

std::ostream& operator<<(std::ostream& out_, const CmdLine& cmdLine_)
    {
      out_ << cmdLine_.ref_argv[0];
      for (int i=1;i<cmdLine_.ref_argc;++i) 
	{
	  out_ << " " << cmdLine_.ref_argv[i];
	}
      return out_;
    };



#endif
