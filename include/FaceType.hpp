#pragma once
#include <array>
#include <iostream>  
struct FaceType
{
public: typedef enum 
  {
    Triangle = 0,
    Quadrilateral
  } enum_t;   

public: static constexpr unsigned int NumKinds = 2;
public: static constexpr std::array<enum_t,NumKinds> All{{ Triangle, Quadrilateral}};

public: static constexpr const std::array<unsigned int,NumKinds> NumNodes{{3,4}};    
public: static constexpr unsigned int GetNumNodes(const enum_t faceKind_)
  {      
    return NumNodes[faceKind_];
  };

public: template <const enum_t _faceKind> static constexpr unsigned int GetNumNodes()
  {      
    return NumNodes[_faceKind];
  };

};

constexpr std::array<unsigned int,FaceType::NumKinds> FaceType::NumNodes;
constexpr std::array<FaceType::enum_t,2> FaceType::All;
  
inline FaceType::enum_t& operator ++(FaceType::enum_t& self_) noexcept
{
  self_=static_cast<FaceType::enum_t>(self_+1);
  return self_;
};


namespace std
{
  
  inline ostream& operator<<(ostream&s_,
			     const FaceType::enum_t& faceType_) noexcept
  {
    switch(faceType_)
      {
      case FaceType::Triangle:
	{
	  s_ << "Triangle";
	  break;
	}
	
      case FaceType::Quadrilateral:
	{
	  s_ << "Quadrilateral";
	  break;
	}	  
      }
    return s_;
  };
  

}
