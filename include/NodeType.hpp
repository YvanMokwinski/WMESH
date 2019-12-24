#pragma once
  
struct NodeType
{
public: typedef enum 
  {
    Node = 0
  } enum_t;   
    
public: static constexpr std::array<enum_t,1> All{{ Node }};
};

constexpr std::array<NodeType::enum_t,1> NodeType::All;
  
inline NodeType::enum_t& operator ++(NodeType::enum_t& self_) noexcept
{
  self_=static_cast<NodeType::enum_t>(self_+1);
  return self_;
};

namespace std
{
  inline ostream& operator<<(ostream&s_,
			     const NodeType::enum_t& nodeType_) noexcept
  {
    switch(nodeType_)
      {
      case NodeType::Node:
	{
	  s_ << "Node";
	  break;
	}
      }
    return s_;
  };

};
