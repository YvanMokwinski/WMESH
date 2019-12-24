#ifndef __header_Debug_h__
#define __header_Debug_h__

#ifndef NDEBUG
#include "Exceptions/ExceptionDebug.h"
#include "Math.h"
struct Debug
{



  static void IsGreaterOrEqualThan(const string&trace_,const long double & value_,const long double& lowerBound_)
  {
    if ( (value_<lowerBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsGreaterOrEqualThan failed " << std::endl;
	s << "        - from IsGreaterOrEqualThan("<< value_ <<","<< lowerBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected lower bound : " << lowerBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsGreaterOrEqualThan(const string&trace_,const double & value_,const double& lowerBound_)
  {
    if ( (value_<lowerBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsGreaterOrEqualThan failed " << std::endl;
	s << "        - from IsGreaterOrEqualThan("<< value_ <<","<< lowerBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected lower bound : " << lowerBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsGreaterOrEqualThan(const string&trace_,const float & value_,const float& lowerBound_)
  {
    if ( (value_<lowerBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsGreaterOrEqualThan failed " << std::endl;
	s << "        - from IsGreaterOrEqualThan("<< value_ <<","<< lowerBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected lower bound : " << lowerBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };




  static void IsGreaterThan(const string&trace_,const long double & value_,const long double& lowerBound_)
  {
    if ( (value_<=lowerBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsGreaterThan failed " << std::endl;
	s << "        - from IsGreaterThan("<< value_ <<","<< lowerBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected lower bound : " << lowerBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsGreaterThan(const string&trace_,const double & value_,const double& lowerBound_)
  {
    if ( (value_<=lowerBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsGreaterThan failed " << std::endl;
	s << "        - from IsGreaterThan("<< value_ <<","<< lowerBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected lower bound : " << lowerBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsGreaterThan(const string&trace_,const float & value_,const float& lowerBound_)
  {
    if ( (value_<=lowerBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsGreaterThan failed " << std::endl;
	s << "        - from IsGreaterThan("<< value_ <<","<< lowerBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected lower bound : " << lowerBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };








  static void IsLowerOrEqualThan(const string&trace_,const long double & value_,const long double& upperBound_)
  {
    if ( (value_>upperBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsLowerOrEqualThan failed " << std::endl;
	s << "        - from IsLowerOrEqualThan("<< value_ <<","<< upperBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected upper bound : " << upperBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsLowerOrEqualThan(const string&trace_,const double & value_,const double& upperBound_)
  {
    if ( (value_>upperBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsLowerOrEqualThan failed " << std::endl;
	s << "        - from IsLowerOrEqualThan("<< value_ <<","<< upperBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected upper bound : " << upperBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsLowerOrEqualThan(const string&trace_,const float & value_,const float& upperBound_)
  {
    if ( (value_>upperBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsLowerOrEqualThan failed " << std::endl;
	s << "        - from IsLowerOrEqualThan("<< value_ <<","<< upperBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected upper bound : " << upperBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };




  static void IsLowerThan(const string&trace_,const long double & value_,const long double& upperBound_)
  {
    if ( (value_>=upperBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsLowerThan failed " << std::endl;
	s << "        - from IsLowerThan("<< value_ <<","<< upperBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected upper bound : " << upperBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsLowerThan(const string&trace_,const double & value_,const double& upperBound_)
  {
    if ( (value_>=upperBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsLowerThan failed " << std::endl;
	s << "        - from IsLowerThan("<< value_ <<","<< upperBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected upper bound : " << upperBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsLowerThan(const string&trace_,const float & value_,const float& upperBound_)
  {
    if ( (value_>=upperBound_)  )
      {
	std::ostringstream s;
	s << "Debug::IsLowerThan failed " << std::endl;
	s << "        - from IsLowerThan("<< value_ <<","<< upperBound_ <<")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "           - expected upper bound : " << upperBound_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };








  static void IsGreaterThan(const string&trace_,const unsigned int & a_,const unsigned int& b_)
  {
    if ( (a_<=b_)  )
      {
	std::ostringstream s;
	s << "Debug::IsGreaterThan failed " << std::endl;
	s << "        - from IsGreaterThan("<< a_ <<","<< b_ <<")" << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsInRange(const string&trace_,const unsigned int & value_,const unsigned int& lowerBound_,const unsigned int& upperBound_)
  {
    if ( (value_<lowerBound_) || (value_>upperBound_) )
      {
	std::ostringstream s;
	s << "Debug::IsInRange failed " << std::endl;
	s << "        - from IsInRange("<< value_ <<","<< lowerBound_ <<","<< upperBound_ << ")" << std::endl;
	s << "           - expected range : [" << lowerBound_ << "," << upperBound_ << "]" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsInRange(const string&trace_,const char & value_,const char& lowerBound_,const char& upperBound_)
  {
    if ( (value_<lowerBound_) || (value_>upperBound_) )
      {
	std::ostringstream s;
	s << "Debug::IsInRange failed " << std::endl;
	s << "        - from IsInRange("<< value_ <<","<< lowerBound_ <<","<< upperBound_ << ")" << std::endl;
	s << "           - expected range : [" << lowerBound_ << "," << upperBound_ << "]" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsInRange(const string&trace_,const short int & value_,const short int& lowerBound_,const short int& upperBound_)
  {
    if ( (value_<lowerBound_) || (value_>upperBound_) )
      {
	std::ostringstream s;
	s << "Debug::IsInRange failed " << std::endl;
	s << "        - from IsInRange("<< value_ <<","<< lowerBound_ <<","<< upperBound_ << ")" << std::endl;
	s << "           - expected range : [" << lowerBound_ << "," << upperBound_ << "]" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsInRange(const string&trace_,const int & value_,const int& lowerBound_,const int& upperBound_)
  {
    if ( (value_<lowerBound_) || (value_>upperBound_) )
      {
	std::ostringstream s;
	s << "Debug::IsInRange failed " << std::endl;
	s << "        - from IsInRange("<< value_ <<","<< lowerBound_ <<","<< upperBound_ << ")" << std::endl;
	s << "           - expected range : [" << lowerBound_ << "," << upperBound_ << "]" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsInRange(const string&trace_,const long int & value_,const long int& lowerBound_,const long int& upperBound_)
  {
    if ( (value_<lowerBound_) || (value_>upperBound_) )
      {
	std::ostringstream s;
	s << "Debug::IsInRange failed " << std::endl;
	s << "        - from IsInRange("<< value_ <<","<< lowerBound_ <<","<< upperBound_ << ")" << std::endl;
	s << "           - expected range : [" << lowerBound_ << "," << upperBound_ << "]" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsInRange(const string&trace_,const long long int & value_,const long long int& lowerBound_,const long long int& upperBound_)
  {
    if ( (value_<lowerBound_) || (value_>upperBound_) )
      {
	std::ostringstream s;
	s << "Debug::IsInRange failed " << std::endl;
	s << "        - from IsInRange("<< value_ <<","<< lowerBound_ <<","<< upperBound_ << ")" << std::endl;
	s << "           - expected range : [" << lowerBound_ << "," << upperBound_ << "]" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };



  static void IsPositive(const string&trace_,const char & value_)
  {
    if ( value_<= ((char)0) )
      {
	std::ostringstream s;
	s << "Debug::IsPositive failed " << std::endl;
	s << "        - from IsPositive("<< value_ << ")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsPositive(const string&trace_,const short int & value_)
  {
    if ( value_<= ((short int)0) )
      {
	std::ostringstream s;
	s << "Debug::IsPositive failed " << std::endl;
	s << "        - from IsPositive("<< value_ << ")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsPositive(const string&trace_,const int & value_)
  {
    if ( value_<= ((int)0) )
      {
	std::ostringstream s;
	s << "Debug::IsPositive failed " << std::endl;
	s << "        - from IsPositive("<< value_ << ")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsPositive(const string&trace_,const long int & value_)
  {
    if ( value_<= ((long int)0) )
      {
	std::ostringstream s;
	s << "Debug::IsPositive failed " << std::endl;
	s << "        - from IsPositive("<< value_ << ")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsPositive(const string&trace_,const long long int & value_)
  {
    if ( value_<= ((long long int)0) )
      {
	std::ostringstream s;
	s << "Debug::IsPositive failed " << std::endl;
	s << "        - from IsPositive("<< value_ << ")" << std::endl;
	s << "           - value : " << value_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };



  static void AreEqual(const string&trace_,const double& a_,const double& b_,const double tolerance = Math<double>::machineEpsilon)
  {
    const R x = (a_>b_) ? a_-b_ : b_-a_;
    if (x > tolerance)
      {
	std::ostringstream s;
	s << "Debug::AreEqual failed " << std::endl;
	s << "        - from AreEqual("<< a_ <<","<< b_ << ")" << std::endl;
	s << "            - expected value : " << a_ << std::endl;
	s << "            - value          : " << b_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };


  static void AreEqual(const string&trace_,const unsigned int& a_,const unsigned int& b_)
  {
    if (a_!=b_)
      {
	std::ostringstream s;
	s << "Debug::AreEqual failed " << std::endl;
	s << "        - from AreEqual("<< a_ <<","<< b_ << ")" << std::endl;
	s << "            - expected value : " << a_ << std::endl;
	s << "            - value          : " << b_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };


  static void AreEqual(const string&trace_,const float& a_,const float& b_,const float tolerance = Math<float>::machineEpsilon)
  {
    const R x = (a_>b_) ? a_-b_ : b_-a_;
    if (x > tolerance)
      {
	std::ostringstream s;
	s << "Debug::AreEqual failed " << std::endl;
	s << "        - from AreEqual("<< a_ <<","<< b_ << ")" << std::endl;
	s << "            - expected value : " << a_ << std::endl;
	s << "            - value          : " << b_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void AreEqual(const string&trace_,const bool& a_,const bool& b_)
  {
    if (a_!=b_)
      {
	std::ostringstream s;
	s << "Debug::AreEqual failed " << std::endl;
	s << "        - from AreEqual("<< a_ <<","<< b_ << ")" << std::endl;
	s << "            - expected value : " << a_ << std::endl;
	s << "            - value          : " << b_ << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsTrue(const string&trace_,const bool& a_)
  {
    if (!a_)
      {
	std::ostringstream s;
	s << "Debug::IsTrue failed " << std::endl;
	s << "        - from IsTrue("<< a_ <<")" << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

  static void IsNotNull(const string&trace_,const void * a_)
  {
    if (!a_)
      {
	std::ostringstream s;
	s << "Debug::IsNotNull failed " << std::endl;
	s << "        - from IsNotNull("<< a_ <<")" << std::endl;
	s << "        - trace :" << trace_;
	throw(ExceptionDebug(s.str()));
      }
  };

};
#endif

#endif
