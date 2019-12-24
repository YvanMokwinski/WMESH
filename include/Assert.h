#ifndef __header_Assert_h__
#define __header_Assert_h__
#include "Exceptions/ExceptionAssert.h"
#include "Math.h"

struct Assert
{
  static void ExpectedAnException(const string&trace_=string("no specified trace"))
  {
    throw(Exception( string("Assert::ExpectedAnException failed, trace : ") + trace_));
  };


  static void IsNotNull(const void * a,const string&trace_=string("no specified trace"))
  {
    if (!a)
      {
	std::cerr << "Assert::IsNotNull(const void *) failed" << std::endl;
	throw(ExceptionAssert(string("Assert::IsNotNull failed, trace: ") + trace_));
      }
  };

  static void IsTrue(const bool& a,const string&trace_=string("no specified trace"))
  {
    if (!a)
      {
	std::cerr << "Assert::IsTrue(bool) failed" << std::endl;
	throw(ExceptionAssert(string("Assert::IsTrue failed, trace: ") + trace_));
      }
  };

  static void IsFalse(const bool& a,const string&trace_=string("no specified trace"))
  {
    if (a)
      {
	std::cerr << "Assert::IsFalse(bool) failed" << std::endl;
	throw(ExceptionAssert(string("Assert::IsFalse failed, trace: ") + trace_));
      }
  };

  static void IsTrue(const size_t& a,const string&trace_=string("no specified trace"))
  {
    if (0!=a)
      {
	std::cerr << "Assert::IsTrue(size_t) failed" << std::endl;
	throw(ExceptionAssert(string("Assert::IsTrue failed, trace: ") + trace_));
      }
  };

  static void AreEqual(const void* a,const void* b,const string&trace_=string("no specified trace"))
  {
    if (a!=b)
      {
	std::cerr << "Assert::AreEqual(const void*,const void*) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(const void*,const void*) failed, trace: ") + trace_));
      }
  };

  static void AreEqual(const bool& a,const bool& b,const string&trace_=string("no specified trace"))
  {
    if (a!=b)
      {
	std::cerr << "Assert::AreEqual(bool,bool) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(bool,bool) failed, trace: ") + trace_));
      }
  };

  static void AreEqual(const string& a,const string& b,const string&trace_=string("no specified trace"))
  {
    if (a!=b)
      {
	std::cerr << "Assert::AreEqual(string,string) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(string,string) failed, trace: ") + trace_));
      }
  };

  static void AreEqual(const unsigned int& a,const unsigned int& b,const string&trace_=string("no specified trace"))
  {
    if (a!=b)
      {
	std::cerr << "Assert::AreEqual(unsigned int,unsigned int) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(unsigned int,unsigned int) failed, trace: ") + trace_));
      }
  };

  static void AreEqual(const int& a,const int& b,const string&trace_=string("no specified trace"))
  {
    if (a!=b)
      {
	std::cerr << "Assert::AreEqual(int, int) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(int, int) failed, trace: ") + trace_));
      }
  };


  static void AreEqual(const long int& a,const long int& b,const string&trace_=string("no specified trace"))
  {
    if (a!=b)
      {
	std::cerr << "Assert::AreEqual(long int, long int) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(long int, long int) failed, trace: ") + trace_));
      }
  };


  static void AreEqual(const long long int& a,const long long int& b,const string&trace_=string("no specified trace"))
  {
    if (a!=b)
      {
	std::cerr << "Assert::AreEqual(long long int, long long int) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(long long int, long long int) failed, trace: ") + trace_));
      }
  };



  static void AreEqual(const long double& a,const long double& b,const string&trace_=string("no specified trace"),const long double tolerance = Math<long double>::machineEpsilon)
  {
    const long double x = (a>b) ? a-b : b-a;
    if (x > tolerance)
      {
	std::cerr << "Assert::AreEqual(long double,long double) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	std::cerr << "        diff is        : " << Math<long double>::Abs(a-b) << std::endl;
	std::cerr << "        tolerance is   : " << tolerance << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(long double,long double) failed, trace: ") + trace_));
      }
  };

  static void AreEqual(const double& a,const double& b,const string&trace_=string("no specified trace"),const double tolerance  = Math<double>::machineEpsilon)
  {
    const double x = (a>b) ? a-b : b-a;
    if (x > tolerance)
      {
	std::cerr << "Assert::AreEqual(double,double) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	std::cerr << "        diff is        : " << Math<double>::Abs(a-b) << std::endl;
	std::cerr << "        tolerance is   : " << tolerance << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(double,double) failed, trace: ") + trace_));
      }
  };

  static void AreEqual(const float& a,const float& b,const string&trace_=string("no specified trace"),const float tolerance = Math<float>::machineEpsilon)
  {
    const float x = (a>b) ? a-b : b-a;
    if (x > tolerance)
      {
	std::cerr << "Assert::AreEqual(float,float) failed" << std::endl;
	std::cerr << "        expected value : " << a << std::endl;
	std::cerr << "        but was        : " << b << std::endl;
	std::cerr << "        diff is        : " << Math<float>::Abs(a-b) << std::endl;
	std::cerr << "        tolerance is   : " << tolerance << std::endl;
	throw(ExceptionAssert(string("Assert::AreEqual(float,float) failed, trace: ") + trace_));
      }
  };
};


#endif
