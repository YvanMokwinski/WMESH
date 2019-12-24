#pragma once
#include <array>

struct VolumeType
{
public: typedef enum 
  {
    Tetrahedron = 0,
    Pyramid,
    Wedge,
    Hexahedron
  } enum_t;
  
public: static constexpr unsigned int NumKinds = 4;
  
public: static constexpr std::array<enum_t,NumKinds> All{{Tetrahedron,Pyramid,Wedge,Hexahedron}};
  
public: static constexpr const std::array<unsigned int,NumKinds> NumNodes{{4,5,6,8}};
  
public: static constexpr unsigned int GetNumNodes(const enum_t volumeKind_)
  {      
    return NumNodes[volumeKind_];
  };

#if 0  
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
#endif    

  
  static constexpr unsigned int s_edgesToNodesTetrahedron[6*2] = { 1,2,
								   2,0,
								   0,1,
								   2,3,
								   0,3,
								   1,3};
  
  static constexpr unsigned int s_edgesToNodesHexahedron[12*2] = { 0,1,
								   1,2,
								   2,3,
								   3,0,
								   4,5,
								   5,6,
								   6,7,
								   7,4,
								   0,4,
								   1,5,
								   2,6,
								   3,7};
  
  static constexpr unsigned int s_edgesToNodesWedge[9*2] = { 0,1,
							     1,2,
							     2,0,
							     3,4,
							     4,5,
							     5,3,
							     0,3,
							     1,4,
							     2,5};
  
  static constexpr unsigned int s_edgesToNodesPyramid[8*2] = { 0,1,
							       1,2,
							       2,3,
							       3,0,
							       0,4,
							       1,4,
							       2,4,
							       3,4};
  
  static constexpr const unsigned int * EdgesToNodes[4] = {s_edgesToNodesTetrahedron,
							   s_edgesToNodesPyramid,
							   s_edgesToNodesWedge,
							   s_edgesToNodesHexahedron};

  static constexpr const unsigned int* GetEdgeToNodes(const enum_t 		volumeKind_,
						      const unsigned int 	localEdgeIndex_)
  {      
    return &EdgesToNodes[volumeKind_][2*localEdgeIndex_];
  };


  static constexpr unsigned int s_facesToNodesTetrahedron[4*4] 
  = { 1,2,3,0,
      2,0,3,0,
      0,1,3,0,
      0,2,1,0};

  static constexpr unsigned int s_facesToNodesHexahedron[4*6] 
  = { 0,3,2,1,
      4,5,6,7,	  
      0,1,5,4,
      1,2,6,5,
      2,3,7,6,
      3,0,4,7};
  
  static constexpr unsigned int s_facesToNodesPyramid[4*5] = { 0,3,2, 1,
							       0,1,4, 0,
							       1,2,4, 0,
							       2,3,4, 0,
							       3,0,4, 0};
  
  static constexpr unsigned int s_facesToNodesWedge[4*5] = { 0,1,2, 0,
							     3,5,4, 0,
							     
							     0,2,5,3,
							     4,5,2,1,
							     0,3,4,1};
  
  static constexpr const unsigned int * FacesToNodes[4] = {s_facesToNodesTetrahedron,
							   s_facesToNodesPyramid,
							   s_facesToNodesWedge,
							   s_facesToNodesHexahedron};
  
  static constexpr const unsigned int* GetFaceToNodes(const enum_t 		volumeKind_,
						      const unsigned int 	localFaceIndex_)
  { 
    return &FacesToNodes[volumeKind_][4*localFaceIndex_];
  };
 
  
  static constexpr const FaceType::enum_t s_faceTypesTetrahedron[4] = {FaceType::Triangle,
								       FaceType::Triangle,
								       FaceType::Triangle,
								       FaceType::Triangle};
  static constexpr const FaceType::enum_t s_faceTypesPyramid[5] = {FaceType::Quadrilateral,
								   FaceType::Triangle,
								   FaceType::Triangle,
								   FaceType::Triangle,
								   FaceType::Triangle};
  static constexpr const FaceType::enum_t s_faceTypesWedge[5] = {  FaceType::Triangle,
								   FaceType::Triangle,
								   FaceType::Quadrilateral,
								   FaceType::Quadrilateral,
								   FaceType::Quadrilateral};
  
  static constexpr const FaceType::enum_t s_faceTypesHexahedron[6] = {  FaceType::Quadrilateral,
									FaceType::Quadrilateral,
									FaceType::Quadrilateral,
									FaceType::Quadrilateral,
									FaceType::Quadrilateral,
									FaceType::Quadrilateral};

  static constexpr const FaceType::enum_t  * FaceTypes[4] = {s_faceTypesTetrahedron,
							     s_faceTypesPyramid,
							     s_faceTypesWedge,
							     s_faceTypesHexahedron};
  
  static constexpr const FaceType::enum_t  GetFaceType(const enum_t 	volumeKind_,
						       const unsigned int 	localFaceIndex_)
  { 
    return FaceTypes[volumeKind_][localFaceIndex_];
  };

  
  
};

constexpr std::array<unsigned int, VolumeType::NumKinds> VolumeType::NumNodes;
constexpr std::array<VolumeType::enum_t, VolumeType::NumKinds> VolumeType::All;
  
inline VolumeType::enum_t & operator ++(VolumeType::enum_t& self_) noexcept { self_ = static_cast<VolumeType::enum_t>(self_+1); return self_; };

#if 0
namespace std
{
  inline ostream& operator<<(ostream&s_,
			     const VolumeType::enum_t& volumeType_) noexcept
  {
#if 0
    switch(volumeType_)
      {
	
      case VolumeType::Tetrahedron:
	{
	  s_ << "Tetrahedron";
	  break;
	}
	  
      case VolumeType::Pyramid:
	{
	  s_ << "Pyramid";
	  break;
	}
	  
      case VolumeType::Wedge:
	{
	  s_ << "Wedge";
	  break;
	}
	  
      case VolumeType::Hexahedron:
	{
	  s_ << "Hexahedron";
	  break;
	}
	  
      }
#endif
    return s_;
  };
};

#endif
