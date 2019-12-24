#pragma once
#include "ReferenceShapeVolume.hpp"
#include "ReferenceShapeFaceTriangle.hpp"
#include "ReferenceShapeFaceQuadrilateral.hpp"

template <>  struct ReferenceShape<DimensionType::Volume,VolumeType::Pyramid>
{ 
public: static constexpr const DimensionType::enum_t Dimension = DimensionType::Volume;
public: static constexpr const unsigned int NbNodes = 5;
public: static constexpr const unsigned int NbEdges = 8;
public: static constexpr const unsigned int NbFaces = 5;
public: static constexpr const unsigned int NbTriangleFaces = 4;
public: static constexpr const unsigned int NbQuadrilateralFaces = 1;
    
public: static constexpr const std::array<unsigned int,Dimension+1> NumEntities{{NbNodes,NbEdges,NbFaces,1}};

public: static constexpr const std::array<FaceType::enum_t, NbFaces>
FaceTypes
  {{ FaceType::Quadrilateral,
	FaceType::Triangle,
	FaceType::Triangle,
	FaceType::Triangle,
	FaceType::Triangle}};
    
public: static constexpr const std::array<unsigned int,NbTriangleFaces> TriangleLocalFaceIndices
  { {1,2,3,4} };
    
public: static constexpr const std::array< unsigned int,NbQuadrilateralFaces> QuadrilateralLocalFaceIndices
  { {0} };
    
public: static constexpr const std::array< std::array<unsigned int, 3> ,NbTriangleFaces> TrianglesToNodes
  {
    std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{0,1,4}},
      std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{1,2,4}},
	std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{2,3,4}},
	  std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{3,0,4}}
  };
    
public: static constexpr const std::array< std::array<unsigned int,4> ,NbQuadrilateralFaces> QuadrilateralsToNodes
  {
    std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{0,3,2,1}}
  };
    
    
public: static constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes> , NbEdges> EdgesToNodes
  {
    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,1}},
      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,2}},
	std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,3}},
	  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{3,0}},
	    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,4}},
	      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,4}},
		std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,4}},
		  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{3,4}}
  };
    
};

using ReferenceShapeVolumePyramid	= ReferenceShapeVolume<VolumeType::Pyramid>;
  

constexpr const std::array<unsigned int,ReferenceShapeVolumePyramid::Dimension+1>
  ReferenceShapeVolumePyramid::NumEntities;  
  
constexpr const std::array<FaceType::enum_t,ReferenceShapeVolumePyramid::NbFaces>
ReferenceShapeVolumePyramid::FaceTypes;  
  
constexpr const std::array<unsigned int, ReferenceShapeVolumePyramid::NbTriangleFaces>
ReferenceShapeVolumePyramid::TriangleLocalFaceIndices;
  
constexpr const std::array< std::array<unsigned int, ReferenceShapeFaceTriangle::NbNodes > , ReferenceShapeVolumePyramid::NbTriangleFaces >
ReferenceShapeVolumePyramid::TrianglesToNodes;
  
constexpr const std::array< unsigned int, ReferenceShapeVolumePyramid::NbQuadrilateralFaces>
ReferenceShapeVolumePyramid::QuadrilateralLocalFaceIndices;
  
constexpr const std::array< std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes> ,ReferenceShapeVolumePyramid::NbQuadrilateralFaces >
ReferenceShapeVolumePyramid::QuadrilateralsToNodes;
  
constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes> , ReferenceShapeVolumePyramid::NbEdges >
ReferenceShapeVolumePyramid::EdgesToNodes;

