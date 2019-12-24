#include "Point.hpp"
#include "Common/Assert.hpp"
#include <iostream>


template <typename _derivedClass> class UnitTest
{
  //!
  //! @brief Constant static cast of the derived class.
  //!  
private: inline const _derivedClass & asImp() const { return static_cast<const _derivedClass&>(*this); }    
  //!
  //! @brief Static cast of the derived class.
  //!  
private: inline _derivedClass & asImp()  { return static_cast<_derivedClass&>(*this); }    

public: inline void run()
  {
    asImp().run();
  };
  
};


template <long unsigned int _size,typename _real_t>
class PointTests : public UnitTest< PointTests<_size,_real_t> >
{
private:
  using real_t = _real_t;
public:
  
  inline void run() 
  {
    TestGetSize();
    TestSetCod();
    TestGetCod();
    TestRandomAccess();
  };
  
  inline void TestGetSize() const 
  {
    Point<_size,real_t> point;
    Assert::AreEqual(_size,
		     point.size());    
  };

  inline void TestSetCod() const 
  {
    Point<_size,real_t> point;
    point.SetCod(32);    
    Assert::AreEqual((int_t)32,
		     point.GetCod());    
    point.SetCod(7);    
    Assert::AreEqual((int_t)7,
		     point.GetCod());    
  };
  
  inline void TestGetCod() const 
  {
    Point<_size,real_t> point;

    //
    // Check the value by default.
    //
    Assert::AreEqual((int_t)0,
		     point.GetCod());    

    
    point.SetCod(64);    
    Assert::AreEqual((int_t)64,
		     point.GetCod());    
  };

  
  inline void TestRandomAccess() const 
  {    
    Point<_size,real_t> point;

    for (unsigned int i=0;i<_size;++i)
      {
	point[i] = real_t(_size-i);
      }

    for (unsigned int i=0;i<_size;++i)
      {
	Assert::AreEqual(real_t(_size-i), point[i]);
      }
  };

};

template <typename _real_t>
class Point2dTests : public UnitTest<Point2dTests<_real_t> >
{
private:
  using real_t = _real_t;
public:

  inline void run() 
  {
    PointTests<2,real_t> pointTests;
    pointTests.run();
    this->TestConstructor();
  };
  
  inline void TestConstructor() const 
  {
    static constexpr real_t x = real_t(1.0);
    static constexpr real_t y = real_t(1.5);
    Point2d<real_t> point(x, y);
    
    Assert::AreEqual(x, point[0]);
    Assert::AreEqual(y, point[1]);    
  };

};

template <typename _real_t>
class Point3dTests : public UnitTest<Point3dTests<_real_t> >
{
private:
  using real_t = _real_t;
public:
  
  inline void run() 
  {
    PointTests<3,real_t> pointTests;
    pointTests.run();
    this->TestConstructor();
  };
  
  inline void TestConstructor() const 
  {
    static constexpr real_t x = real_t(1.0);
    static constexpr real_t y = real_t(1.5);
    static constexpr real_t z = real_t(3.5);
    Point3d<real_t> point(x, y, z);

    Assert::AreEqual(x, point[0]);
    Assert::AreEqual(x, point[1]);
    Assert::AreEqual(z, point[2]);    
  };

};


class TestFixture
{
private:  unsigned long int numTests = 0;
private:  unsigned long int numTestsFailed = 0;
public:
  ~TestFixture()
  {
    validation();
  };
  
private:  template <typename _unitTestImpl> void RunTest(UnitTest<_unitTestImpl>& unitTest_)
  {
    unitTest_.run();
  };

private: void validation()
  {  
    if (numTestsFailed)
      {
	std::ostringstream s;
	s << "UnitTests failed" << std::endl;
	s << "        NumTests   : " << numTests << std::endl;
	s << "        Failed     : " << numTestsFailed << std::endl;
	throw(Exception(s.str()));     
      }
  };
  
public:  template <typename _unitTest> void test()
  {
    
    ++numTests;
    try
      {
	_unitTest test;
	RunTest(test);
      }
    
    catch(Exception& e)
      {
	std::cerr << e.what() << std::endl;
	++numTestsFailed;
      }
  };
  
};

int main()
{

  
  TestFixture testFixture;

  testFixture.test< PointTests<8,float> >();
  testFixture.test< PointTests<12,float> >();
  testFixture.test< PointTests<17,float> >();

  testFixture.test< PointTests<8,double> >();
  testFixture.test< PointTests<12,double> >();
  testFixture.test< PointTests<17,double> >();

  testFixture.test< PointTests<8,long double> >();
  testFixture.test< PointTests<12,long double> >();
  testFixture.test< PointTests<17,long double> >();

  testFixture.test< Point2dTests<float> >();
  testFixture.test< Point2dTests<double> >();
  testFixture.test< Point2dTests<long double> >();

  testFixture.test< Point3dTests<float> >();
  testFixture.test< Point3dTests<double> >();
  testFixture.test< Point3dTests<long double> >();
  
  return 0;
}
