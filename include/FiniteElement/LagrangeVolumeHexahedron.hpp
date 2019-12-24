#pragma once

#include "LagrangeFaceQuadrilateral.hpp"
#include "LagrangeVolume.hpp"

namespace FiniteElement
{
  
  template <unsigned int _degree> struct Lagrange<_degree,
						  DimensionType::Volume,
						  VolumeType::Hexahedron>
  {
#if 0    
  private: static constexpr unsigned int s_numVerticesInCell 	 = ReferenceShapeVolumeHexahedron::NbNodes;
  private: static constexpr unsigned int s_numEdgesInCell 	 = ReferenceShapeVolumeHexahedron::NbEdges;
  private: static constexpr unsigned int s_numFacesInCell 	 = ReferenceShapeVolumeHexahedron::NbFaces;
    
  public: static constexpr unsigned int NumDofsPerCell 	 	= (_degree+1)*(_degree+1)*(_degree+1);
  public: static constexpr unsigned int NumDofsPerInteriorEdge 	= (_degree>0) ? _degree-1 : 0;
    
  public: static constexpr unsigned int NumDofsPerInteriorFace[2]
    { 0 , LagrangeFaceQuadrilateral<_degree>::NumDofsPerEntities[DimensionType::Face] };
    
  public: static constexpr unsigned int NumDofsPerInteriorVolume = NumDofsPerCell
		   - s_numVerticesInCell
		   - s_numEdgesInCell * NumDofsPerInteriorEdge
		   - 2 * NumDofsPerInteriorFace[0]
		   - 3 * NumDofsPerInteriorFace[1];
#endif

  private: static constexpr unsigned int s_numVerticesInCell 	 = ReferenceShapeVolumeHexahedron::NbNodes;
  private: static constexpr unsigned int s_numEdgesInCell 	 = ReferenceShapeVolumeHexahedron::NbEdges;
  private: static constexpr unsigned int s_numFacesInCell 	 = ReferenceShapeVolumeHexahedron::NbFaces;
  public: static constexpr unsigned int Degree = _degree;
  public: static constexpr VolumeType::enum_t VolumeKind = VolumeType::Hexahedron;
  public: static constexpr unsigned int NumDofs 	 	= (_degree+1)*(_degree+1)*(_degree+1);
  public: static constexpr const unsigned int NumDofsPerInteriorNode[1]{ LagrangeNode<_degree>::NumDofsPerEntities[0]};
  public: static constexpr const unsigned int NumDofsPerInteriorEdge[1]{ LagrangeEdge<_degree>::NumDofsPerEntities[1]};
  public: static constexpr const unsigned int NumDofsPerInteriorFace[2]{ 0, LagrangeFaceQuadrilateral<_degree>::NumDofsPerEntities[2]};
    
  public: static constexpr const unsigned int NumDofsPerInteriorVolume[4]{NumDofs
	- s_numVerticesInCell
	- s_numEdgesInCell * LagrangeEdge<_degree>::NumDofsPerEntities[1]
	- s_numFacesInCell * LagrangeFaceQuadrilateral<_degree>::NumDofsPerEntities[2],0,0,0};
    
  public: static constexpr const unsigned int * NumDofsPerEntities[4]{
    NumDofsPerInteriorNode,
      NumDofsPerInteriorEdge,
      NumDofsPerInteriorFace,
      NumDofsPerInteriorVolume};
    
  };
  
  template <unsigned int _degree> using LagrangeVolumeHexahedron = LagrangeVolume<_degree,
										  VolumeType::Hexahedron>;
  
};  
