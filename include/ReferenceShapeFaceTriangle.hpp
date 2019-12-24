#pragma once
#include "ReferenceShapeFace.hpp"
  
template <>  struct ReferenceShape<DimensionType::Face,FaceType::Triangle>
{
public: static constexpr const DimensionType::enum_t Dimension = DimensionType::Face;
public: static constexpr const unsigned int NbNodes = 3;
public: static constexpr const unsigned int NbEdges = 3;
public: static constexpr const std::array<unsigned int,3> NumEntities{{NbNodes,NbEdges,1}};
public: static constexpr const std::array< std::array<unsigned int, 2 > , 3 >
EdgesToNodes
  {
    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,1}},
      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,2}},
	std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,0}}
  };
};
  
using ReferenceShapeFaceTriangle 		= ReferenceShapeFace<FaceType::Triangle>;
  

constexpr const std::array<unsigned int, ReferenceShapeFaceTriangle::Dimension + 1 > ReferenceShapeFaceTriangle::NumEntities;
constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes> , ReferenceShapeFaceTriangle::NbEdges> ReferenceShapeFaceTriangle::EdgesToNodes;
