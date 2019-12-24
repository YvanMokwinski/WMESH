#pragma once

#include "LagrangeEdge.hpp"

namespace FiniteElement
{
  template <unsigned int _degree,FaceType::enum_t _faceType> using LagrangeFace = Lagrange<_degree,
											   DimensionType::Face,
											   _faceType>;
  
};  
