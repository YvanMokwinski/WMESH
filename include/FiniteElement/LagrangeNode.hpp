#pragma once

#include "Lagrange.hpp"

namespace FiniteElement
{

  template <unsigned int _degree>
  struct Lagrange<_degree,
		  DimensionType::Node,
		  NodeType::Node>
  {
  public: static constexpr unsigned int Degree = _degree;
  public: static constexpr const unsigned int NumDofs = (_degree > 0) ? 1 : 0;
  public: static constexpr const std::array<unsigned int,1> NumDofsPerEntities
    {{ NumDofs }};    
  };

  template <unsigned int _degree> using LagrangeNode = Lagrange<_degree,
								DimensionType::Node,
								NodeType::Node>;
  
};  
