#pragma once
#include "ReferenceShape.hpp"

template <>  struct ReferenceShape<DimensionType::Node,NodeType::Node>
{
public: static constexpr const DimensionType::enum_t Dimension = DimensionType::Node;
public: static constexpr const std::array<unsigned int, Dimension+1> NumEntities{{1}};
};
  
using ReferenceShapeNode = ReferenceShape<DimensionType::Node,NodeType::Node>;
#ifndef NDEBUG
constexpr const std::array<unsigned int,ReferenceShapeNode::Dimension+1> ReferenceShapeNode::NumEntities;
#endif
