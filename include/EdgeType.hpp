#pragma once
  
struct EdgeType
{
public: typedef enum 
  {
    Edge = 0
  } enum_t;   
    
public: static constexpr std::array<enum_t,1> All{{ Edge }};

};

#ifndef NDEBUG
constexpr std::array<EdgeType::enum_t,1> EdgeType::All;
#endif
  
inline EdgeType::enum_t& operator ++(EdgeType::enum_t& self_) noexcept
{
  self_=static_cast<EdgeType::enum_t>(self_+1);
  return self_;
};


namespace std
{
  
  inline ostream& operator<<(ostream&s_,
			     const EdgeType::enum_t& edgeType_) noexcept
  {
    switch(edgeType_)
      {
      case EdgeType::Edge:
	{
	  s_ << "Edge";
	  break;
	}
      }
    return s_;
  };
  
};  
  
