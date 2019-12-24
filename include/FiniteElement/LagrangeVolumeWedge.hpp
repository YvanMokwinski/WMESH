#pragma once

#include "LagrangeFaceTriangle.hpp"
#include "LagrangeFaceQuadrilateral.hpp"
#include "LagrangeVolume.hpp"

namespace FiniteElement
{

  
  template <unsigned int _degree> struct Lagrange<_degree,DimensionType::Volume,VolumeType::Wedge>
  {

  private: static constexpr unsigned int s_numVerticesInCell 	 = ReferenceShapeVolumeWedge::NbNodes;
  private: static constexpr unsigned int s_numEdgesInCell 	 = ReferenceShapeVolumeWedge::NbEdges;
  private: static constexpr unsigned int s_numFacesInCell 	 = ReferenceShapeVolumeWedge::NbFaces;
  
    //
    // 2 * n = (_degree+1)*(_degree+1)*(_degree+1) + (_degree+1)*(_degree+1)
    // n = ( (_degree+1)*(_degree+1)*(_degree+1) + (_degree+1)*(_degree+1) ) / 2
    //
  
  public: static constexpr unsigned int NumDofsPerCell 	 	= ((_degree+1)*(_degree+1)*(_degree+2)) / 2;
  public: static constexpr unsigned int NumDofsPerInteriorEdge 	= (_degree>0) ? _degree-1 : 0;
  
  public: static constexpr unsigned int NumDofsPerInteriorFace[]
    { LagrangeFaceTriangle<_degree>::NumDofsPerEntities[DimensionType::Face],
	LagrangeFaceQuadrilateral<_degree>::NumDofsPerEntities[DimensionType::Face] };
  
  public: static constexpr unsigned int NumDofsPerInteriorVolume = NumDofsPerCell
		   - s_numVerticesInCell
		   - s_numEdgesInCell * NumDofsPerInteriorEdge
		   - 2 * NumDofsPerInteriorFace[0]
		   - 3 * NumDofsPerInteriorFace[1];

  public: static constexpr unsigned int Degree = _degree;
  public: static constexpr VolumeType::enum_t VolumeKind = VolumeType::Wedge;
  public: static constexpr const unsigned int NumDofs = ((_degree+1)*(_degree+1)*(_degree+2)) / 2;
    
  };
  
  template <unsigned int _degree> using LagrangeVolumeWedge = LagrangeVolume<_degree,
									     VolumeType::Wedge>;
  
};  
