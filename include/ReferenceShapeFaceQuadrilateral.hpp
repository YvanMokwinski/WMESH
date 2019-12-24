#pragma once
#include "ReferenceShape.hpp"
  
template <>  struct ReferenceShape<DimensionType::Face,FaceType::Quadrilateral>
{

public: static constexpr const DimensionType::enum_t Dimension = DimensionType::Face;
public: static constexpr const unsigned int NbNodes = 4;
public: static constexpr const unsigned int NbEdges = 4;
public: static constexpr const std::array<unsigned int, Dimension + 1> NumEntities{{NbNodes,NbEdges,1}};
public: static constexpr const std::array< std::array<unsigned int, 2 > , NbEdges >
EdgesToNodes
  {
    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,1}},
      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,2}},
	std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,3}},
	  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{3,0}}
  };
    
};

using ReferenceShapeFaceQuadrilateral	= ReferenceShapeFace<FaceType::Quadrilateral>;

constexpr const std::array<unsigned int, ReferenceShapeFaceQuadrilateral::Dimension + 1 >
  ReferenceShapeFaceQuadrilateral::NumEntities;
  
constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes>, ReferenceShapeFaceQuadrilateral::NbEdges>
ReferenceShapeFaceQuadrilateral::EdgesToNodes;
