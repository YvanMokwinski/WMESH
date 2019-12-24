#pragma once

#include "LagrangeFace.hpp"

namespace FiniteElement
{
  template <unsigned int _degree>
  struct Lagrange<_degree,
		  DimensionType::Face,
		  FaceType::Triangle>
  {    
  public: static constexpr unsigned int Degree = _degree;
  public: static constexpr const unsigned int NumDofs = ((_degree+1)*(_degree+2)) / 2;
  public: static constexpr const std::array<unsigned int,3> NumDofsPerEntities
    {{ LagrangeNode<_degree>::NumDofs,
	  LagrangeEdge<_degree>::NumDofsPerEntities[DimensionType::Edge],
	  (_degree>2) ? ((_degree-1)*(_degree-2))/2 : 0}};
  };  
  template <unsigned int _degree> using LagrangeFaceTriangle = LagrangeFace<_degree,
									    FaceType::Triangle>;
  
};  
