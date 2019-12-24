#pragma once
#include "ReferenceShapeNode.hpp"

template <>  struct ReferenceShape<DimensionType::Edge,EdgeType::Edge>
{
public: static constexpr const DimensionType::enum_t Dimension = DimensionType::Edge;
public: static constexpr const unsigned int NbNodes = 2;
public: static constexpr const std::array<unsigned int, Dimension+1> NumEntities{{NbNodes,1}};
};

using ReferenceShapeEdge = ReferenceShape<DimensionType::Edge,EdgeType::Edge>;
#ifndef NDEBUG
constexpr const std::array<unsigned int,ReferenceShapeEdge::Dimension+1> ReferenceShapeEdge::NumEntities;
#endif
