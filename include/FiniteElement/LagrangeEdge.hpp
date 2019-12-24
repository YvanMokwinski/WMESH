#pragma once

#include "LagrangeNode.hpp"

namespace FiniteElement
{
  
  template <unsigned int _degree>
  struct Lagrange<_degree,
		  DimensionType::Edge,
		  EdgeType::Edge>
  {
  public: static constexpr unsigned int Degree = _degree;
  public: static constexpr const unsigned int NumDofs = _degree+1;
  public: static constexpr const std::array<unsigned int,2> NumDofsPerEntities
    {{ LagrangeNode<_degree>::NumDofs,
	  (_degree>0) ? _degree-1 : 0 }};    
  };

  template <unsigned int _degree> using LagrangeEdge = Lagrange<_degree,
								DimensionType::Edge,
								EdgeType::Edge>;
  
};  
