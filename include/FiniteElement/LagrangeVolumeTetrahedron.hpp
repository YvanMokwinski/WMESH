#pragma once

#include "LagrangeFaceTriangle.hpp"
#include "LagrangeVolume.hpp"

namespace FiniteElement
{

  
  template <unsigned int _degree>
  struct Lagrange<_degree,DimensionType::Volume,VolumeType::Tetrahedron>
  {
#if 0    
  private: static constexpr unsigned int s_numVerticesInCell 	 = ReferenceShapeVolumeTetrahedron::NbNodes;
  private: static constexpr unsigned int s_numEdgesInCell 	 = ReferenceShapeVolumeTetrahedron::NbEdges;
  private: static constexpr unsigned int s_numFacesInCell 	 = ReferenceShapeVolumeTetrahedron::NbFaces;
  
  public: static constexpr unsigned int NumDofsPerCell 	 	= ((_degree+1)*(_degree+2)*(_degree+3)) / 6;
  public: static constexpr unsigned int NumDofsPerInteriorEdge 	= (_degree>0) ? _degree-1 : 0; 
  public: static constexpr unsigned int NumDofsPerInteriorFace[2]
    { LagrangeFaceTriangle<_degree>::NumDofsPerEntities[DimensionType::Face],
	0};
    
  public: static constexpr unsigned int
  NumDofsPerInteriorVolume = NumDofsPerCell
		   - s_numVerticesInCell
		   - s_numEdgesInCell * NumDofsPerInteriorEdge
		   - s_numFacesInCell * NumDofsPerInteriorFace[0];
#endif
    
  private: static constexpr unsigned int s_numVerticesInCell 	 = ReferenceShapeVolumeTetrahedron::NbNodes;
  private: static constexpr unsigned int s_numEdgesInCell 	 = ReferenceShapeVolumeTetrahedron::NbEdges;
  private: static constexpr unsigned int s_numFacesInCell 	 = ReferenceShapeVolumeTetrahedron::NbFaces;

  public: static constexpr unsigned int Degree = _degree;
  public: static constexpr VolumeType::enum_t VolumeKind = VolumeType::Tetrahedron;
  public: static constexpr unsigned int NumDofs 	 	= ((_degree+1)*(_degree+2)*(_degree+3)) / 6;  
  public: static constexpr const unsigned int NumDofsPerInteriorNode[1]{ LagrangeNode<_degree>::NumDofsPerEntities[0]};
  public: static constexpr const unsigned int NumDofsPerInteriorEdge[1]{ LagrangeEdge<_degree>::NumDofsPerEntities[1]};
  public: static constexpr const unsigned int NumDofsPerInteriorFace[2]{ LagrangeFaceTriangle<_degree>::NumDofsPerEntities[2],
	0};
    
  public: static constexpr const unsigned int NumDofsPerInteriorVolume[4]{NumDofs
	- s_numVerticesInCell
	- s_numEdgesInCell * LagrangeEdge<_degree>::NumDofsPerEntities[1]
	- s_numFacesInCell * LagrangeFaceTriangle<_degree>::NumDofsPerEntities[2],0,0,0};
  public: static constexpr const unsigned int * NumDofsPerEntities[4]{
    NumDofsPerInteriorNode,
      NumDofsPerInteriorEdge,
      NumDofsPerInteriorFace,
      NumDofsPerInteriorVolume};

  };
  
  
  template <unsigned int _degree> using LagrangeVolumeTetrahedron = LagrangeVolume<_degree,
										   VolumeType::Tetrahedron>;
  
};  
