#pragma once
#include "ReferenceShapeVolume.hpp"
#include "ReferenceShapeFaceTriangle.hpp"
#include "ReferenceShapeFaceQuadrilateral.hpp"
  
template <>  struct ReferenceShape<DimensionType::Volume,VolumeType::Wedge>
{ 
public: static constexpr const DimensionType::enum_t Dimension = DimensionType::Volume;
public: static constexpr const unsigned int NbNodes = 6;
public: static constexpr const unsigned int NbEdges = 9;
public: static constexpr const unsigned int NbFaces = 5;
public: static constexpr const unsigned int NbTriangleFaces = 2;
public: static constexpr const unsigned int NbQuadrilateralFaces = 3;
    
public: static constexpr const std::array<unsigned int,Dimension+1> NumEntities{{NbNodes,NbEdges,NbFaces,1}};
public: static constexpr const std::array<FaceType::enum_t, NbFaces>
FaceTypes
  {{ FaceType::Triangle,
	FaceType::Triangle,
	FaceType::Quadrilateral,
	FaceType::Quadrilateral,
	FaceType::Quadrilateral}};
    
public: static constexpr const std::array<unsigned int,NbTriangleFaces> TriangleLocalFaceIndices
  { {0,1} };
    
public: static constexpr const std::array< unsigned int,NbQuadrilateralFaces> QuadrilateralLocalFaceIndices
  { {2,3,4} };
    
public: static constexpr const std::array< std::array<unsigned int, ReferenceShapeFaceTriangle::NbNodes> ,NbTriangleFaces> TrianglesToNodes
  {
    std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{0,1,2}},
      std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{3,5,4}}
  };
    
public: static constexpr const std::array< std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes> ,NbQuadrilateralFaces> QuadrilateralsToNodes
  {
    std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{0,2,5,3}},
      std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{4,5,2,1}},
	std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes>{{0,3,4,1}}
  };

    
public: static constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes> , NbEdges> EdgesToNodes
  {
    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,1}},
      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,2}},
	std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,0}},
	  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{3,4}},
	    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{4,5}},
	      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{5,3}},
		std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,3}},
		  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,4}},
		    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,5}}
  };
    
    
};

using ReferenceShapeVolumeWedge 		= ReferenceShapeVolume<VolumeType::Wedge>;
  


constexpr const std::array<unsigned int,ReferenceShapeVolumeWedge::Dimension+1>
  ReferenceShapeVolumeWedge::NumEntities;  
  
constexpr const std::array<FaceType::enum_t,ReferenceShapeVolumeWedge::NbFaces>
ReferenceShapeVolumeWedge::FaceTypes;  
  
constexpr const std::array<unsigned int, ReferenceShapeVolumeWedge::NbTriangleFaces>
ReferenceShapeVolumeWedge::TriangleLocalFaceIndices;
  
constexpr const std::array< std::array<unsigned int, ReferenceShapeFaceTriangle::NbNodes > , ReferenceShapeVolumeWedge::NbTriangleFaces >
ReferenceShapeVolumeWedge::TrianglesToNodes;
  
constexpr const std::array< unsigned int, ReferenceShapeVolumeWedge::NbQuadrilateralFaces>
ReferenceShapeVolumeWedge::QuadrilateralLocalFaceIndices;
  
constexpr const std::array< std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes> ,ReferenceShapeVolumeWedge::NbQuadrilateralFaces >
ReferenceShapeVolumeWedge::QuadrilateralsToNodes;

constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes> , ReferenceShapeVolumeWedge::NbEdges >
ReferenceShapeVolumeWedge::EdgesToNodes;

  
