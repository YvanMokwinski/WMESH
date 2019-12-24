#include <iostream>
#include "GenericEncoding.hpp"
#include "Common/Assert.hpp"

template <typename _int_t,int _nbits> class GenericEncodingTests
{
  using this_t = GenericEncoding<_int_t,_nbits>;
  
public: template <typename T = _int_t>
static inline typename std::enable_if<std::is_signed<T>::value>::type Run()
  {
    
    //
    // testing from -100 to 100
    //
    for (T up=-100;up<=100;++up)
      {
	for (T low = 0;low < 4;++low)
	  {
	    auto encodedValue = this_t::Encod(up, low);	    
	    Assert::AreEqual(low,this_t::Low(encodedValue));
	    Assert::AreEqual(up,this_t::Up(encodedValue));
	  }
      }
    
    // std::cout << "testing ["<< this_t::s_upLimit - 3 << "," << this_t::s_upLimit << "[" << std::endl;    
    // std::cout << "testing from "<< this_t::s_upLimit - 3 << " to " << this_t::s_upLimit << std::endl;    
    for (T up=this_t::s_upLimit-3;up < this_t::s_upLimit;++up)
      {
	for (T low = 0;low < 4;++low)
	  {
	    auto encodedValue = this_t::Encod(up, low);
	    Assert::AreEqual(low,this_t::Low(encodedValue));
	    Assert::AreEqual(up,this_t::Up(encodedValue));	    
	  }
      }

    // std::cout << "testing ["<< -(this_t::s_upLimit - 1) << "," << -(this_t::s_upLimit-10) << "[" << std::endl;    
    for (T up=-(this_t::s_upLimit-1);up < -(this_t::s_upLimit-10);++up)
      {	
	for (T low = 0;low < 4;++low)
	  {
	    auto encodedValue = this_t::Encod(up,low);
	    Assert::AreEqual(low,this_t::Low(encodedValue));
	    Assert::AreEqual(up,this_t::Up(encodedValue));	    
	  }
      }
  };

  
public: template <typename T = _int_t>
static inline typename std::enable_if<!std::is_signed<T>::value >::type Run()
  {    

    // std::cout << "testing from 0 to 100"  << std::endl;
    for (T up=0;up<=100;++up)
      {
	for (T low = 0;low < 4;++low)
	  {
	    const auto encodedValue = this_t::Encod(up, low);
	    Assert::AreEqual(low,this_t::Low(encodedValue));
	    Assert::AreEqual(up,this_t::Up(encodedValue));
	    if (this_t::Low(encodedValue) != low)
	      {
		std::cout << "low is wrong, expected " << low << " but got " <<  this_t::Low(encodedValue) << std::endl;
		exit(1);
	      }
	    if (this_t::Up(encodedValue) != up)
	      {
		std::cout << "up is wrong, expected " << up << " but got " <<  this_t::Up(encodedValue) << std::endl;
		exit(1);
	      }
	    
	  }
      }

    // std::cout << "testing ["<< this_t::s_upLimit - 3 << "," << this_t::s_upLimit << "[" << std::endl;    
    // std::cout << "testing from "<< this_t::s_upLimit - 3 << " to " << this_t::s_upLimit << std::endl;    
    for (T up=this_t::s_upLimit-3;up < this_t::s_upLimit;++up)
      {
	for (T low = 0;low < 4;++low)
	  {
	    const auto encodedValue = this_t::Encod(up, low);
	    Assert::AreEqual(low,this_t::Low(encodedValue));
	    Assert::AreEqual(up,this_t::Up(encodedValue));
	  }
      }

  };
};

int main()
{
  GenericEncodingTests<int,2>::Run();
  GenericEncodingTests<unsigned int,2>::Run();

  GenericEncodingTests<long int,2>::Run();
  GenericEncodingTests<unsigned long int,2>::Run();

  GenericEncodingTests<short int,2>::Run();
  GenericEncodingTests<unsigned short int,2>::Run();

  return 0;
}
