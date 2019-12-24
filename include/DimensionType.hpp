#pragma once

struct DimensionType
{

public: typedef enum
  { Node=0,
    Edge=1,
    Face=2,
    Volume=3 }
    enum_t;
  
public:
  static constexpr std::array<enum_t,4> All{{Node,Edge,Face,Volume}};
};

#ifndef NDEBUG
constexpr std::array<DimensionType::enum_t,4> DimensionType::All;
#endif
  
DimensionType::enum_t & operator ++(DimensionType::enum_t& self_)
{
  self_=static_cast<DimensionType::enum_t>(self_+1);
  return self_;
};

namespace std
{
  std::ostream& operator<<(std::ostream&s_,
			   const DimensionType::enum_t& d_)
  {
    switch(d_)
      {
      case DimensionType::Node:
	{
	  s_ << "Node";
	  break;
	}

      case DimensionType::Edge:
	{
	  s_ << "Edge";
	  break;
	}

      case DimensionType::Face:
	{
	  s_ << "Face";
	  break;
	}

      case DimensionType::Volume:
	{
	  s_ << "Volume";
	  break;
	}
	  
      }
    return s_;
  };
  
};
