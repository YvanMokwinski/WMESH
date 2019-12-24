#pragma once
#include "ReferenceShapeVolume.hpp"
#include "ReferenceShapeFaceTriangle.hpp"
  
template <>  struct ReferenceShape< DimensionType::Volume, VolumeType::Tetrahedron >
{
public: static constexpr const DimensionType::enum_t Dimension = DimensionType::Volume;
public: static constexpr const unsigned int NbNodes = 4;
public: static constexpr const unsigned int NbEdges = 6;
public: static constexpr const unsigned int NbFaces = 4;
public: static constexpr const unsigned int NbTriangleFaces = 4;
public: static constexpr const unsigned int NbQuadrilateralFaces = 0;
    
public: static constexpr const std::array<unsigned int,Dimension + 1> NumEntities{{NbNodes,NbEdges,NbFaces,1}};
    
public: static constexpr const std::array<FaceType::enum_t, NbTriangleFaces > FaceTypes
  {
    {FaceType::Triangle,
	FaceType::Triangle,
	FaceType::Triangle,
	FaceType::Triangle}
  };

public: static constexpr const std::array<unsigned int, NbTriangleFaces >
TriangleLocalFaceIndices { {0,1,2,3} };
    
public: static constexpr const std::array< unsigned int, NbQuadrilateralFaces >
QuadrilateralLocalFaceIndices { };
    
public: static constexpr const std::array< std::array<unsigned int, ReferenceShapeFaceTriangle::NbNodes > , NbTriangleFaces >
TrianglesToNodes
  {
    std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{1,2,3}},
      std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{2,0,3}},
	std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{0,1,3}},
	  std::array<unsigned int,ReferenceShapeFaceTriangle::NbNodes>{{0,2,1}}
  };
    
public: static constexpr const std::array< std::array<unsigned int, NbTriangleFaces > , NbQuadrilateralFaces >
QuadrilateralsToNodes{};
    
public: static constexpr const std::array< std::array<unsigned int, 2 > , NbEdges >
EdgesToNodes
  {
    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,2}},
      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,0}},
	std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,1}},
	  std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{2,3}},
	    std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{0,3}},
	      std::array<unsigned int,ReferenceShapeEdge::NbNodes>{{1,3}}
  };
    
};


using ReferenceShapeVolumeTetrahedron 	= ReferenceShapeVolume<VolumeType::Tetrahedron>;
  

constexpr const std::array<unsigned int,ReferenceShapeVolumeTetrahedron::Dimension+1>
  ReferenceShapeVolumeTetrahedron::NumEntities;  
  
constexpr const std::array<FaceType::enum_t,ReferenceShapeVolumeTetrahedron::NbFaces>
ReferenceShapeVolumeTetrahedron::FaceTypes;  
  
constexpr const std::array<unsigned int, ReferenceShapeVolumeTetrahedron::NbTriangleFaces>
ReferenceShapeVolumeTetrahedron::TriangleLocalFaceIndices;
  
constexpr const std::array< std::array<unsigned int, ReferenceShapeFaceTriangle::NbNodes > , ReferenceShapeVolumeTetrahedron::NbTriangleFaces >
ReferenceShapeVolumeTetrahedron::TrianglesToNodes;
  
constexpr const std::array< unsigned int, ReferenceShapeVolumeTetrahedron::NbQuadrilateralFaces>
ReferenceShapeVolumeTetrahedron::QuadrilateralLocalFaceIndices;
  
constexpr const std::array< std::array<unsigned int,ReferenceShapeFaceQuadrilateral::NbNodes> ,ReferenceShapeVolumeTetrahedron::NbQuadrilateralFaces >
ReferenceShapeVolumeTetrahedron::QuadrilateralsToNodes;

constexpr const std::array< std::array<unsigned int,ReferenceShapeEdge::NbNodes> , ReferenceShapeVolumeTetrahedron::NbEdges >
ReferenceShapeVolumeTetrahedron::EdgesToNodes;

