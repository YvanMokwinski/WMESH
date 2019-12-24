#pragma once

#include "LagrangeFace.hpp"

namespace FiniteElement
{
  
  template <unsigned int _degree>
  struct Lagrange<_degree,
		  DimensionType::Face,
		  FaceType::Quadrilateral>
  {    
  public: static constexpr unsigned int Degree = _degree;
  public: static constexpr const unsigned int NumDofs = (_degree+1)*(_degree+1);
  public: static constexpr const std::array<unsigned int,3> NumDofsPerEntities
    {{ LagrangeNode<_degree>::NumDofs,
	  LagrangeEdge<_degree>::NumDofsPerEntities[DimensionType::Edge],
	  (_degree-1)*(_degree-1)}};
  };  
  template <unsigned int _degree> using LagrangeFaceQuadrilateral = LagrangeFace<_degree,
										 FaceType::Quadrilateral>;
  
};  
