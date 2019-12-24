#pragma once
#include "ReferenceShapeVolume.hpp"
#include "ReferenceShapeFaceQuadrilateral.hpp"

template <>  struct ReferenceShape<DimensionType::Volume,VolumeType::Hexahedron>
{
public: static constexpr const DimensionType::enum_t Dimension = DimensionType::Volume;
public: static constexpr const unsigned int NbNodes = 8;
public: static constexpr const unsigned int NbEdges = 12;
public: static constexpr const unsigned int NbFaces = 6;
public: static constexpr const unsigned int NbTriangleFaces = 0;
public: static constexpr const unsigned int NbQuadrilateralFaces = NbFaces;
    
public: static constexpr const std::array<unsigned int, Dimension + 1> NumEntities{{NbNodes,NbEdges,NbFaces,1}};
    
public: static constexpr const std::array<FaceType::enum_t,NbFaces>  FaceTypes
  {{ FaceType::Quadrilateral,
	FaceType::Quadrilateral,
	FaceType::Quadrilateral,
	FaceType::Quadrilateral,
	FaceType::Quadrilateral,
	FaceType::Quadrilateral}};
    
public: static constexpr const std::array<unsigned int, NbTriangleFaces>
TriangleLocalFaceIndices { };
    
public: static constexpr const std::array< std::array<unsigned int, ReferenceShapeFaceTriangle::NbNodes > , 0>
TrianglesToNodes { };
    
public: static constexpr const std::array< unsigned int, NbQuadrilateralFaces >
QuadrilateralLocalFaceIndices { {0,1,2,3,4,5} };
    
public: static constexpr const std::array< std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes> , NbFaces >
QuadrilateralsToNodes
  {
    std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{0,3,2,1}},
      std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{4,5,6,7}},	  
	std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{0,1,5,4}},
	  std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{1,2,6,5}},
	    std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{2,3,7,6}},
	      std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{3,0,4,7}}
  };
    
public: static constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes> , NbEdges >
EdgesToNodes
  {
    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,1}},
      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,2}},
	std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,3}},
	  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{3,0}},
	    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{4,5}},
	      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{5,6}},
		std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{6,7}},
		  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{7,4}},
		    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,4}},
		      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,5}},
			std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,6}},
			  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{3,7}}
  }; 
    
};

using ReferenceShapeVolumeHexahedron 		= ReferenceShapeVolume<VolumeType::Hexahedron>;
  

constexpr const std::array<unsigned int,ReferenceShapeVolumeHexahedron::Dimension+1>
  ReferenceShapeVolumeHexahedron::NumEntities;  
  
constexpr const std::array<FaceType::enum_t,ReferenceShapeVolumeHexahedron::NbFaces>
ReferenceShapeVolumeHexahedron::FaceTypes;  
  
constexpr const std::array<unsigned int, ReferenceShapeVolumeHexahedron::NbTriangleFaces>
ReferenceShapeVolumeHexahedron::TriangleLocalFaceIndices;
  
constexpr const std::array< std::array<unsigned int, ReferenceShapeFaceTriangle::NbNodes > , ReferenceShapeVolumeHexahedron::NbTriangleFaces >
ReferenceShapeVolumeHexahedron::TrianglesToNodes;
  
constexpr const std::array< unsigned int, ReferenceShapeVolumeHexahedron::NbQuadrilateralFaces>
ReferenceShapeVolumeHexahedron::QuadrilateralLocalFaceIndices;
  
constexpr const std::array< std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes> ,ReferenceShapeVolumeHexahedron::NbQuadrilateralFaces >
ReferenceShapeVolumeHexahedron::QuadrilateralsToNodes;

constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes> , ReferenceShapeVolumeHexahedron::NbEdges >
ReferenceShapeVolumeHexahedron::EdgesToNodes;
